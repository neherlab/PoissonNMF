import ij.*;
import ij.process.*;
import ij.gui.*;
import ij.io.*;

import java.awt.*;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.util.Arrays;
import java.io.*;
import java.util.StringTokenizer;
import java.util.Vector;

import ij.plugin.filter.*;

/**
 * PoissonNMF_ 
 * imagej plugin for poissonNMF
 *
 * @author <a href="mailto:fabian.theis@helmholtz-muenchen.de">Fabian Theis</a>
 * @author <a href="mailto:neher@kitp.ucsb.edu">Richard Neher</a>
 */
public class PoissonNMF_ implements PlugInFilter {
	String PoissonNMFversion="0.8.8";
	
	ImagePlus imp;	//input, holds the stack X
	ImageStack X; 	//emission signal
	ImagePlus img;	//output, hold the stack A
	ImageStack A;	//concentrations 
	ImagePlus modifiedimg;	//output, hold the stack A
	ImageStack modifiedA;	//concentrations 
	
	//dimension of the images
	int r=3;
	int n;
	int w; //width;
	int h; //height;
	int wh; //product
	int tz_dim; //number of product of dimensions 4 and 5
	int [] dim=new int [5];

	//internal error codes
	int CANCELLED=-1;
	int NOROISELECTED=-2;
	int GOODINPUT=1;
	int BADVALUE=-3;
	
	//matrices that hold the NMF objects and intermediates
	int []bgMask;
	double[][] Xsub;	//signal
	double[][] Asub;	//concentration
	double[][] ASsub;	//expected concentration
	double[][] S;		//spectra
	double[][] modifiedS;		//spectra
	float[][] initialS;	//initial spectra
	double[][] pinvS;	//pseudo inverse of spectra matrix
	
	//auxillary quantities and parameters
	double[][] channel_lambdas;	//center wavelength of the channels
	int [] channel_order;
	int [] inverse_channel_order;
	double channel_width;		//channel width
	double[][] Segbiasterm;		//temporary variable to store gradient contribution of the seg bias
	double[] bg; 		//background signal (vector, one entry for each spectral channel)
	double[] bg_sigma;	//std of the above
	double[] power; 	//buffer
	float signal_nothing=1;		//positive but irrelevant signal strength
	double conc_nothing=1;		//positive but irrelevant concentration
	double E=1; 		//E=1 one-norm for spectra, E=2 two norm
	double segbias=0;	//weight of the segregation bias
	int noPix;			//number of pixels in the current subsample
	boolean isHyperstack;
	double bg_threshold=50;		// Threshold, below which data is not used  
	double no_std=2;				// number of standard deviations of background noise, below which data is not used
	double saturation_threshold=4000;	//threshold above which the data is discarded
	boolean []spectra_fixed;			//vector of booleans to decide which spectra to fix during optimization

	int maxit=100;						//number of iterations
	int subsamples=3;					//subsample levels, each with 10fold less data
	int no_preiterations=10;				// number of iterations on concentrations without spectra update
	int no_postiterations=5;				// number of iterations on concentrations without spectra update before showing concentrations
	boolean cancelled=false;
	boolean interrupted=false;
	String [] spec_choice;			//array of strings with the choice of initial spectra
	String bg_choice;				//choice for the background mode
	PlotWindow plotw=null;			//window handle to display the spectra
	String datatype;
	String [] datatype_choices;
	String config_file;
	/**
	 * inital stuff
	 * @param arg
	 * @param imp
	 * @return
	 * @see ij.plugin.filter.PlugInFilter#setup(java.lang.String, ij.ImagePlus)
	 */
	public int setup(String arg, ImagePlus imp_in) {
		imp=imp_in;
		datatype_choices=new String [3]; 
		datatype_choices[0]="Regular stack";
		datatype_choices[1]="Leica SP2";
		datatype_choices[2]="Zeiss LSM";
		datatype=datatype_choices[0];
		if (imp==null)
		{
			GenericDialog datasource= new GenericDialog("Select Data type");
			datasource.addChoice("Type:", datatype_choices, datatype);
			datasource.showDialog();
			if (datasource.wasCanceled())
			{	
				return DONE;
			}
			else
			{
				datatype=datasource.getNextChoice();
				if (datatype.equals(datatype_choices[1]))
					{if (openLeicaSP2()==CANCELLED) return DONE;}
				else
					{openStack();}
			}
		}
		if (imp!=null)
		{
			//get image dimensions and check whether there are several image slices
			dim=imp.getDimensions();
			int d=1; for (int i = 2; i < dim.length; i++) d*=dim[i];
			if (d==1)
			{
				IJ.error("Stack required!");
				return DONE; 	
			}
			//if stack is ok, convert to 32 gray scale and assign dimensions.
			if ((imp.getBytesPerPixel()<32)) new StackConverter(imp).convertToGray32();
			X = imp.getStack();
			isHyperstack=imp.isHyperStack();
			w = dim[0]; //X.getWidth();
			h = dim[1]; //X.getHeight();
			wh = w*h;	//number of pixels per frame
						//for hyperstacks, save number of frames in tz_dim
			if (isHyperstack) {n = dim[2]; tz_dim=dim[3]*dim[4];}
			else {n=X.getSize(); tz_dim=1;}
			//allocate memory to specify which points are above the background threshold
			bgMask=new int [wh*tz_dim];
		}
		else
		{
			IJ.error("Stack required!");
			return DONE; 
		}
		if (arg.equals("about")) { 
			showAbout(); 
			return DONE; 
		} 
		return NO_IMAGE_REQUIRED+DOES_32 + DOES_8G + DOES_16; //just convert the last ones
	}

	/**
	 * main method for running the algorithm
	 * @param ip
	 * @see ij.plugin.filter.PlugInFilter#run(ij.process.ImageProcessor)
	 */
	public void run(ImageProcessor ip) {
		//check whether the version of ImageJ is new enough
		if (IJ.versionLessThan("1.39t"))
		{
			IJ.error("PoissonNMF requires ImageJ version 1.39t or greater!");
			return;
		}
		//user input of basic parameters
		if (parameterDialog()<0) return;

		int []iterations=new int [subsamples];
		double sum=0;
		/**
		 * normalize the segregation bias such that the contribution of the additional term to the cost
		 * function is of the order of that of a background pixel when segbias was 1
		 * the necessary magnitude of the segregation bias  depends on the overlap of the 
		 * spectra and the distribution of dye concentrations
		 */
		segbias*=n*0.5*Math.log(Math.max(bg_threshold,100));
		//calculate the number of iterations made for each subsample level
		for (int s=subsamples-1; s>=0; s--) sum+=Math.pow(2,s);
		for (int s=subsamples-1; s>=0; s--) iterations[s]=(int) (maxit*Math.pow(2,s)/sum);
		double sum2=0;
		for (int s=subsamples-1; s>=0; s--) sum2+=Math.pow(0.1, s)*iterations[s];
		
		//the iteration loop -- the number of iterations is specified by the user
		double count=0;
		PoissonNMFprogress progress=new PoissonNMFprogress(this);

		subsample(Math.pow(0.1, subsamples-1));	//produce a new subsample of size 10^-s
		for (int s=subsamples-1; s>=0; s--)
		{
			subsample(Math.pow(0.1, s));	//produce a new subsample of size 10^-s
			if (s==subsamples-1) 
			{		
				initA();								//initialize the concentrations with random numbers
				for (int it=0; it<no_preiterations; it++) updateA();	//do a small number of concentration updates
			}
			else solveforA();								//initialize the concentrations with least squares solution
			for (int it=0; it<iterations[s]; it++) {
				IJ.showProgress(count/sum2);
				normalize();	//normalize spectra and concentrations
				updateS();		//update the spectra
				updateA();		//update the concentrations
				count+=Math.pow(0.1, s);
				if (cancelled)	return;	//check state of cancel button
				else if (interrupted) { s=-1;  break;}	//break out of loops if interrupted
				plotSpectra();
			}
		}
		if (!(cancelled || interrupted)) progress.closeWindow();
		IJ.showProgress(1);

		//display results
		PoissonNMFspectracorrection speccorr; 
		showConcentrations();
		results_panel();
	}

	
	/**
	 * Function that reads the previous settings from the imageJ preferences file
	 **/
	private void getPrefs() {
		String label;
		maxit=(int)Prefs.get("PoissonNMF.maxit", maxit);
		segbias=Prefs.get("PoissonNMF.segbias", segbias);
		saturation_threshold=Prefs.get("PoissonNMF.saturation_threshold", saturation_threshold);
		bg_threshold=Prefs.get("PoissonNMF.bg_threshold", bg_threshold);
		bg_choice=Prefs.get("PoissonNMF.bg_choice", "none");
		subsamples=(int)Prefs.get("PoissonNMF.subsamples", subsamples);

		//choice of initial spectra and decision to keep some spectra fixed
		if ((int)Prefs.get("PoissonNMF.r",0)==r)
		{
			for(int dye=0; dye<r; dye++)
			{			
				label="PoissonNMF.Dye_";
				label=label.concat(Integer.toString(dye+1));
				spec_choice[dye]=Prefs.get(label, "none");
				label="PoissonNMF.DyeFixed_";
				label=label.concat(Integer.toString(dye+1));
				spectra_fixed[dye]=Prefs.get(label, false);
			}
		}
		//channel boundaries
		double lambda_max=650;
		double lambda_min=480;
		channel_width=(lambda_max-lambda_min)/(n-1.0);
		channel_lambdas[0][0]=lambda_min-0.5*channel_width;
		channel_lambdas[0][1]=lambda_min+0.5*channel_width;
		//set them to the default values
		for (int i = 1; i < channel_lambdas.length; i++) {
			channel_lambdas[i][0]=channel_lambdas[i-1][0]+channel_width;
			channel_lambdas[i][1]=channel_lambdas[i-1][1]+channel_width;
		}
		//Override default values if previous values are available
		if ((int)Prefs.get("PoissonNMF.n",0)==n)
		{
			for (int i = 0; i < channel_lambdas.length; i++) {
				label="PoissonNMF.Channel_lower_";
				label=label.concat(Integer.toString(i+1));
				channel_lambdas[i][0]=Prefs.get(label, channel_lambdas[i][0]);
				label="PoissonNMF.Channel_upper_";
				label=label.concat(Integer.toString(i+1));
				channel_lambdas[i][1]=Prefs.get(label, channel_lambdas[i][1]);
			}
		}

	}
	/**
	 * Save the parameters to the imageJ preferences file: one-to-one correspondence
	 * to getPrefs() above
	 */
	private void setPrefs() {
		String label;
		Prefs.set("PoissonNMF.maxit", maxit);
		Prefs.set("PoissonNMF.segbias", segbias);
		Prefs.set("PoissonNMF.saturation_threshold", saturation_threshold);
		Prefs.set("PoissonNMF.bg_threshold", bg_threshold);
		Prefs.set("PoissonNMF.bg_choice", bg_choice);
		Prefs.set("PoissonNMF.subsamples", subsamples);
		Prefs.set("PoissonNMF.r",r);
		//initial spectra
		for(int dye=0; dye<r; dye++)
		{			
			label="PoissonNMF.Dye_";
			label=label.concat(Integer.toString(dye+1));
			Prefs.set(label, spec_choice[dye]);
			label="PoissonNMF.DyeFixed_";
			label=label.concat(Integer.toString(dye+1));
			Prefs.set(label, spectra_fixed[dye]);
		}
		//channel boundaries
		Prefs.set("PoissonNMF.n",n);
		for (int i = 0; i < channel_lambdas.length; i++) {
			label="PoissonNMF.Channel_lower_";
			label=label.concat(Integer.toString(i+1));
			Prefs.set(label, channel_lambdas[i][0]);
			label="PoissonNMF.Channel_upper_";
			label=label.concat(Integer.toString(i+1));
			Prefs.set(label, channel_lambdas[i][1]);
		}
	}

	/**
	 * once the dimensions of the problem are known (i.e. the number of 
	 * sources), allocate the necessary matrices
	 */
	private void allocate_spectra()
	{
		S = new double [r][n];
		initialS = new float [r][n];
		spectra_fixed=new boolean [r];		
		power=new double [r];
		bg=new double[n];
		bg_sigma=new double [n];
		spec_choice=new String [r];
		channel_lambdas=new double [n][2];
		channel_order=new int [n];
		inverse_channel_order=new int [n];
	}
	
	
	/**
	 * init S with Gaussians
	 */
	private void gauss_spectra() {
		int y,ii;int z;
		//spectra - Gauss shape
		for (z=0;z<r;z++) {		//loop over dyes
			spectra_fixed[z]=false;
			float sigma = (n-1f)/(r+1f)+.1f;	//variance of gauss spectra
			for (y=1;y<=n;y++) 				//loop over channels and assign spectra
			{
				ii=channel_order[y-1];
				S[z][ii]= (float) (1/(sigma*Math.sqrt(2*Math.PI))*
						Math.exp(-.5*Math.pow(y-1-((z+1)*n)/(r+1f),2)/Math.pow(sigma,2)))*
						(channel_lambdas[ii][1]-channel_lambdas[ii][0]);
			}
		}
	}

	/**
	 * init Asub with random number
	 */
	private void initA() {
		int x;int z;

		for (x=0;x<noPix;x++)		//loop over pixels and assign concentrations with random numbers
			for (z=0;z<r;z++) {		//loop over dyes
				Asub[x][z]= (float) (Math.random()/2+.5);
			}
	}

	/**
	 * init Asub with least square solution using current spectra
	 */
	private void solveforA() {
		int x;int lambda, source;
		calc_pinv(false);
		for (x=0;x<noPix;x++)		//loop over pixels and assign concentrations with random numbers
		{
			for (source=0; source<r; source++)
			{
				Asub[x][source]=0;
				for (lambda=0; lambda<n; lambda++)
				{	
					Asub[x][source]+=pinvS[source][lambda]*Xsub[x][lambda];
				}
			}
		}
	}


	/**
	 * normalize everything (Asub and S) to unit power
	 */
	private void normalize() {
		int x;
		int y;
		int z;
		for (z=0;z<r;z++) {				//for each dye, calculate the 'E'-norm of the spectrum
			for (y=0,power[z]=0;y<n;y++)
				power[z]+=Math.pow(S[z][y],E);
			power[z]=(float) Math.pow(power[z], 1f/E);
		}
		for (z=0;z<r;z++) {			//rescale S and A accordingly
			for (y=0;y<n;y++)
				S[z][y]/=power[z];
			for (x=0;x<noPix;x++)
				Asub[x][z]*=power[z];
		}
	}

	/**
	 * calculate the maxima of each dye in a certain slice of the image  
	 * This is needed to display the RGB overlay, where each source is normalized
	 * to one.
	 * @param currentslice_conc is the slice (time, depth) currently on display
	 */
	private void get_ConcMaxima(float[] maxima, int currentslice_conc) {
		if (maxima.length!=r)
		{
			IJ.error("Bad array length!");
			return;
		}
		int x;
		int z;
		float[] rowA;
		for (z=0;z<r;z++) {
			maxima[z]=0;
			rowA = (float[])A.getPixels(currentslice_conc+z+1);
			for (x=0;x<wh;x++)			//determine max of each image
				if (maxima[z]<rowA[x]) maxima[z]=rowA[x];
		}
	}
	
	/**
	 * choose subsample from above background data points, subsample fraction is frac
	 */
	private void subsample(double frac)
	{
		noPix=0;
		int count=0, currentslice_signal, currentpix;
		//for each pixel with bgMask>0 (not background or saturated), draw a random number
		//if selected, set bgMask to 2n
		//if bgMask is already 2n, nothing happens -> pixels from previous subsamples are always included
		for (int slice=0; slice<tz_dim; slice++) //loop over all t and z slices
		{
			currentpix=slice*wh;
			for (int pix=currentpix; pix<currentpix+wh; pix++)
			{	
				if (bgMask[pix]>0 && Math.random()<frac) {bgMask[pix]=2*n;}
				if (bgMask[pix]==2*n) noPix++;
			}
		}
		Xsub=new double [noPix][n];
		Asub=new double[noPix][r];
		Segbiasterm=new double [noPix][r];
		ASsub=new double[noPix][n];
		//loop over slices and pixels and copy the signal into from the ImagePlus
		//into Xsub, pixels are chosen on the basis of bgMask, count is the continuous index
		for (int slice=0; slice<tz_dim; slice++)
		{
			currentslice_signal=slice*n;
			currentpix=slice*wh;
			for (int pix=0; pix<wh; pix++)
			{
				if (bgMask[currentpix+pix]==2*n)	//if pixel is selected, add it to the data set
				{
					for (int lambda=0; lambda<n; lambda++)
					{	//subtract the background signal, set negative pixels to signal_nothing
						Xsub[count][lambda]= ((float [])X.getPixels(currentslice_signal+lambda+1))[pix]-bg[lambda];
						if (Xsub[count][lambda]<signal_nothing) Xsub[count][lambda]=signal_nothing;
					}
					count++;
				}
			}
		}
	}



	/**
	 * plot/show image stack of concentrations (A)
	 * 
	 */
	private void showConcentrations() {
		int lambda, source, count=0;

		segbias=0;
		for(int i=0; i<no_postiterations; i++) updateA();
		
		//allocate a new image stack, one slice for each dye at each t and z slice
		A=new ImageStack(w,h);
		for (int j=0; j<tz_dim; j++)	//time and depth slices
			for (int i=0; i<r; i++)		//add slice for each dye
				A.addSlice("source "+(i+1), new FloatProcessor(w,h) );


		calc_pinv(false);
		double [] c=new double [r];
		int currentslice_conc, currentslice_signal, currentpix;
		//loop over all slices and pixels and put together the image stack of the NMF sources
		for (int slice=0; slice<tz_dim; slice++)
		{
			currentpix=slice*wh;
			currentslice_conc=slice*r;
			currentslice_signal=slice*n;
			for (int pix=0; pix<wh; pix++)
			{
				if (bgMask[currentpix+pix]==2*n) 	//fill in NMF pixels
				{
					for (source=0; source<r; source++)
						((float[])A.getPixels(currentslice_conc+source+1))[pix]=(float)Asub[count][source];
					count++;
				}
				else	//solve least squares for all other pixels
				{
					//construct vector, set negative values to signal_nothing
					//solve for concentrations
					for (source=0; source<r; source++)
					{
						c[source]=0;
						for (lambda=0; lambda<n; lambda++)
						{	
							c[source]+=pinvS[source][lambda]*Math.max(((float[])X.getPixels(currentslice_signal+lambda+1))[pix]-bg[lambda],signal_nothing);
						}
					}
					for (source=0; source<r; source++)
						((float[])A.getPixels(currentslice_conc+source+1))[pix]=(float) c[source];
				}
			}
		}
		//Display
		img = new ImagePlus("NMF sources", A);
		if (isHyperstack) //for hyperstacks, set the time and depth dimensions
		{
			img.setDimensions(r,dim[3], dim[4]);
			img.setOpenAsHyperStack(true);
		}
		img.show();
		img.updateAndDraw();
	}
	

	/**
	 * plot/show image stack of concentrations (A)
	 * 
	 */
	public void showModifiedConcentrations() {
		int lambda, source;

		segbias=0;
		for(int i=0; i<no_postiterations; i++) updateA();
		
		//allocate a new image stack, one slice for each dye at each t and z slice
		if (modifiedA==null){
			modifiedA=new ImageStack(w,h);
			for (int j=0; j<tz_dim; j++)	//time and depth slices
				for (int i=0; i<r; i++)		//add slice for each dye
					modifiedA.addSlice("source "+(i+1), new FloatProcessor(w,h) );
		}

		calc_pinv(true);
		double [] c=new double [r];
		int currentslice_conc, currentslice_signal;
		//loop over all slices and pixels and put together the image stack of the NMF sources
		for (int slice=0; slice<tz_dim; slice++)
		{
			currentslice_conc=slice*r;
			currentslice_signal=slice*n;
			for (int pix=0; pix<wh; pix++)
			{
				//construct vector, set negative values to signal_nothing
				//solve for concentrations
				for (source=0; source<r; source++)
				{
					c[source]=0;
					for (lambda=0; lambda<n; lambda++)
					{	
						c[source]+=pinvS[source][lambda]*Math.max(((float[])X.getPixels(currentslice_signal+lambda+1))[pix]-bg[lambda],signal_nothing);
					}
				}
				for (source=0; source<r; source++)
					((float[])modifiedA.getPixels(currentslice_conc+source+1))[pix]=(float) c[source];
			}
		}
		//Display
		if (modifiedimg==null){
			modifiedimg = new ImagePlus("Corrected NMF sources", modifiedA);
			if (isHyperstack) //for hyperstacks, set the time and depth dimensions
			{
				modifiedimg.setDimensions(r,dim[3], dim[4]);
				modifiedimg.setOpenAsHyperStack(true);
			}
			modifiedimg.show();
		}
		modifiedimg.updateAndDraw();
	}

	/**
	 * show the RGB overlay of specified sources
	 */
	public void RGB_overlay()
	{
		float [] conc_maxima=new float [r];
		int [] sources=new int [3];
		int currentslice_conc;
		//prompt the user for an assignment of channels to Red, Green, Blue
		GenericDialog sourcesDialog=new GenericDialog("Select Sources");
		if (r>2)
		{
			sourcesDialog.addNumericField("Blue:", 1, 0);
			sourcesDialog.addNumericField("Green:", 2, 0);
			sourcesDialog.addNumericField("Red:", 3, 0);
		}else
		{	
			sourcesDialog.addNumericField("Green:", 1, 0);
			sourcesDialog.addNumericField("Red:", 2, 0);
		}				
		sourcesDialog.showDialog();
		if (sourcesDialog.wasCanceled())
		{	
			return;
		}
		else
		{
			if (isHyperstack)	//if hyperstack, determine the tz-slice on display
			{
				currentslice_conc=img.getCurrentSlice();
				currentslice_conc--;	//correct for base 1 counting
				currentslice_conc-=currentslice_conc%r;	//determine the first channel
			}
			else currentslice_conc=0;
			sources[0]=(int) sourcesDialog.getNextNumber();
			sources[1]=(int) sourcesDialog.getNextNumber();
			if (r>2) sources[2]=(int) sourcesDialog.getNextNumber();
			//check that entered source values are positive and not larger than the number of sources
			for (int i=0; i<Math.min(r,3);i++)
				if (sources[i]>r  || sources[i]<1) 
				{
					IJ.showMessage("Invalid source number!");
					return;
				}
			ImagePlus overlay=null;
			overlay=NewImage.createRGBImage("Overlay of concentrations", w,h,1, NewImage.FILL_BLACK);
			int source;
			//normalize images to max==1
			get_ConcMaxima(conc_maxima, currentslice_conc);
			//produce the RGB overlay if the number of sources == 3
			if (r>2)
			{
				for (int pix=0; pix<wh; pix++)
				{
					for (source=0; source<3; source++)
						((int[])overlay.getProcessor().getPixels())[pix]+=
							((int)(((1<<8)-1.0)*Math.abs(((float[])A.getPixels(currentslice_conc+sources[source]))[pix])/conc_maxima[sources[source]-1]))<<(8*source);
				}
			}
			else
			{
				for (int pix=0; pix<wh; pix++)
				{
					for (source=0; source<2; source++)
						((int[])overlay.getProcessor().getPixels())[pix]+=
							((int)(((1<<8)-1.0)*Math.abs(((float[])A.getPixels(currentslice_conc+sources[source]))[pix])/conc_maxima[sources[source]-1]))<<(8*(source+1));
				}
			}
			overlay.show();
			overlay.updateAndDraw();
		}
	}

	/**
	 * Fucntion that draws a map of the parts of the image that are background, signal and
	 * saturated. meant as a sanity check for the user.
	 */
	public void BG_map()
	{
		int currentpix;
		ImagePlus regions=null;
		regions=NewImage.createRGBImage("Background and saturated regions", w,h,1, NewImage.FILL_BLACK); 
		if (isHyperstack)	//determine the current tz-slice on display
		{
			currentpix=img.getCurrentSlice();
			currentpix--;				//correct for base 1 counting
			currentpix-=currentpix%r;	//determine first color channel slice
			currentpix/=r;				//number of current tz-slice
			currentpix*=wh;				//pixel offset = #tz-slices * pixel per frame
		}
		else currentpix=0;
		//construct the map of background, signal and saturated regions
		for (int pix=0; pix<wh; pix++)
		{
			if (bgMask[currentpix+pix]==2*n) 	//forground pixels
			{
				((int[])regions.getProcessor().getPixels())[pix]=1<<15;
			}
			else	
			{
				if (bgMask[currentpix+pix]>=0) //background pixel
					((int[])regions.getProcessor().getPixels())[pix]=1<<7;
				else							//saturated pixels
					((int[])regions.getProcessor().getPixels())[pix]=1<<23;
			}
		}
		regions.show();
		regions.updateAndDraw();
	}
	
	/**
	 * generate new plot for the spectra
	 * @return 
	 * 
	 */
	public void plotSpectra() {
		int y,ii;
		int z;
		double [][] spec=new double [r][n];
		for (z=0;z<r;z++) {
			for (y=0;y<n;y++) {
				ii=channel_order[y];
				spec[z][y]=S[z][ii];
			}
		}
		//colors for the different spectra
		Color[] colors = {Color.blue,Color.green, Color.red, Color.cyan, Color.gray, Color.darkGray};
		Plot plot = new Plot("PoissonNMF spectra","wave length [nm]","intensity",(float[])null,(float[])null, PlotWindow.LINE);
	
		if (plotw==null){
			plotw=new PlotWindow("PoissonNMF spectra","wave length [nm]","intensity",(float[])null,(float[])null);
		}
		float[] rowS = new float[n];
		float[] lambdas=new float [n];
		for (int i = 0; i < lambdas.length; i++) {
			ii=channel_order[i];
			lambdas[i]=0.5f*((float)(channel_lambdas[ii][0]+channel_lambdas[ii][1])); //channel wavelength as x-axis
		}
		plot.setLimits(channel_lambdas[channel_order[0]][0], channel_lambdas[channel_order[n-1]][1], 0, 1);
		int c=0; float maxS=0;
		//determine maximum of all spectra
		for (z=0;z<r;z++) {
			for (y=0;y<n;y++) {
				ii=inverse_channel_order[y];
				if (spec[z][y]/(channel_lambdas[ii][1]-channel_lambdas[ii][0])>maxS) 
					maxS=(float) ((float) spec[z][y]/(channel_lambdas[ii][1]-channel_lambdas[ii][0]));
			}
		}
		//loop over spectra and plot each
		for (z=0;z<r;z++) {
			for (y=0;y<n;y++) {
				ii=inverse_channel_order[y];
				rowS[y]=(float) (spec[z][y]/(channel_lambdas[ii][1]-channel_lambdas[ii][0])/maxS);
			}
			plot.setLineWidth(2);            
			plot.setColor(colors[c%colors.length]);
			plot.addPoints(lambdas, rowS, PlotWindow.LINE);
			//plot initial spectra only if not kept fixed
			if (spectra_fixed[z]==false)
			{
				plot.setLineWidth(1);            
				plot.addPoints(lambdas, initialS[z], PlotWindow.LINE);
			}
			if (modifiedS!=null)
			{
				for (y=0;y<n;y++) {
					ii=inverse_channel_order[y];
					rowS[y]=(float) (modifiedS[z][y]/(channel_lambdas[ii][1]-channel_lambdas[ii][0])/maxS);
				}
				plot.setLineWidth(3);            
				plot.setColor(colors[c%colors.length]);
				plot.addPoints(lambdas, rowS, PlotWindow.LINE);				
			}
			c++;
		}
		plot.draw();
		plotw.drawPlot(plot);
	}

	/**
	 * Plot background spectrum
	 * @return
	 */
	public PlotWindow plotbg() {
		int y;
		float [] bgplot=new float [n];
		//color fot the spectrum
		Color[] colors = {Color.blue};
		Plot plot=new Plot("Background spectrum","wave length [nm]","intensity",(float[])null,(float[])null,PlotWindow.LINE);
		float[] lambdas=new float [n];
		for (int i = 0; i < lambdas.length; i++) {
			lambdas[i]=0.5f*((float)(channel_lambdas[i][0]+channel_lambdas[i][1])); //channel wavelength as x-axis
		}
		double minbg=saturation_threshold;
		double maxbg=0;
		for (y = 0; y < n; y++) {
			bgplot[y]=(float) bg[y];
			if (bg[y]>maxbg) maxbg=bg[y]; 
			if (bg[y]<minbg) minbg=bg[y]; 
		}
		if (minbg==maxbg) {minbg-=5; maxbg+=5;}
		plot.setLimits(channel_lambdas[0][0], channel_lambdas[n-1][1], minbg, maxbg);

		//loop over spectra and plot each
		plot.setLineWidth(2);            
		plot.setColor(colors[0]);
		plot.addPoints(lambdas, bgplot, PlotWindow.LINE);
		return plot.show();
	}

	/**
	 * Function that loops over all dyes and presents the user a save dialog
	 * to save each sprectrum
	 * @return
	 */
	public int saveSpectra(){
		//Dialog asking which spectra to save
		GenericDialog save_spectra=new GenericDialog("Save some spectra?");
		String [] labels=new String [r];
		Boolean [] spectra_save=new Boolean [r];
		for (int z = 0; z < r; z++) {
			labels[z]= "Save spectrum "+new Integer(z+1).toString()+"?";
			save_spectra.addCheckbox(labels[z], true);
		}
		save_spectra.showDialog();
		if (save_spectra.wasCanceled())
		{	
			return CANCELLED;
		}
		else
		{	//save result in boolean vector
			for (int z = 0; z < r; z++) {
				spectra_save[z]=save_spectra.getNextBoolean();
			}
		}


		for (int i = 0; i < spectra_save.length; i++) {
			if (spectra_save[i])	//save requested spectra
			{
				SaveDialog sd=new SaveDialog("Save spectrum as ...", "spectrum"+new Integer(i+1).toString(), ".emn");
				String dir=sd.getDirectory();
				String file=sd.getFileName();
				if (file==null) return -1;
				try{
					int channel;
					BufferedWriter out=new BufferedWriter(new FileWriter(dir+file));
					for (int l = 0; l < n; l++){	//write wavelength in column 1, spectrum in column 2
						channel=channel_order[l];
						out.write(new Float((channel_lambdas[channel][0]+channel_lambdas[channel][1])*0.5).toString()+"\t"
								+new Float(S[i][channel]/(channel_lambdas[channel][1]-channel_lambdas[channel][0])).toString()+"\t");
						out.write("\n");
					}
					out.close();
				}
				catch (Exception e){
					IJ.error("Error while saving spectra:", e.getMessage());
				}
			}
		}
		return 0;
	}

	public void show_simplex()
	{
		int [] dyes=new int [3];
		if (r<3)
		{
			IJ.showMessage("Simplex plots require at least three dyes!");
			return;
		}
		else if (r==3)
		{
			dyes[0]=0; dyes[1]=1; dyes[2]=2;
 		}
		else if (r>3)
		{
			//make user input the dyes for the projection.
			GenericDialog sourcesDialog=new GenericDialog("Select dyes for simplex projection");
			if (r>2)
			{
				sourcesDialog.addNumericField("Dye 1:", 1, 0);
				sourcesDialog.addNumericField("Dye 2:", 2, 0);
				sourcesDialog.addNumericField("Dye 3:", 3, 0);
			}
			sourcesDialog.showDialog();
			if (sourcesDialog.wasCanceled())
			{	
				return;
			}
			else
			{
				dyes[0]=(int) sourcesDialog.getNextNumber()-1;
				dyes[1]=(int) sourcesDialog.getNextNumber()-1;
				dyes[2]=(int) sourcesDialog.getNextNumber()-1;
				//check whether input is distinct and a valid dye
				if (!(dyes[0]>=0 && dyes[1]>=0 &&dyes[2]>=0 && dyes[0]<r && dyes[1]<r 
						&& dyes[2]<r && dyes[0]!=dyes[1] &&  dyes[0]!=dyes[2] && dyes[1]!=dyes[2]))
				{
					IJ.showMessage("Bad input!");
					return;
				}
			}
		}
		modifiedS=new double [r][n];
		for (int dye=0; dye<r; dye++)
		{
			for (int ch=0; ch<n; ch++){
				modifiedS[dye][ch]=S[dye][ch];				
			}
		}
		PoissonNMFspectracorrection speccorr=new PoissonNMFspectracorrection(this, dyes);
	}
	/**
	 * about dialog
	 */
	void showAbout() { 
		IJ.showMessage("About PoissonNMF...",
				"An ImageJ plugin for separating multiple dyes from an image stack. Version "+PoissonNMFversion
		); 
	} 

	/**
	 * User interface to input background signal, either by ROI or by number 
	 */
	private int getBackground(String choice)
	{	
		if (choice==null)
		{
			IJ.error("Null pointer for background choice!");
			return BADVALUE;
		}
		
		String[] bgitems=new String [4];
		bgitems[0]="Mininal values";
		bgitems[1]="ROI selection";
		bgitems[2]="manually";
		bgitems[3]="flat";

		int lambda;
		int goodinput=0;
		if (choice.equals(bgitems[1]))		//open a dialog asking the user to select a ROI
		{
			new WaitForUserDialog("Background Selection", "Select Background ROI").show();
			Roi	bg_roi=imp.getRoi();

			if (bg_roi==null)
				goodinput=NOROISELECTED;	//no ROI selected -- have the user try again
			else
				goodinput=GOODINPUT;

			if (goodinput==GOODINPUT)	//read ROI, calculate background intensity
			{
				//determine current slice of the image to extract the background
				int currentslice_signal=imp.getCurrentSlice();
				currentslice_signal--;	//correct for counting with base 1
				currentslice_signal-=currentslice_signal%n;	//determine first channel
				int no_bgpix=0; 
				float v=0;
				for (int x=0; x<w; x++) //loop over the entire image (this is terribly inefficient)
				{
					for (int y=0; y<h; y++)
					{
						if (bg_roi.contains(x,y))	//if pixel in ROI, use it to calculate 
							//background mean and sigma for each channel
						{	
							no_bgpix++;
							for (lambda=0; lambda<n; lambda++)
							{
								v=((float [])X.getPixels(currentslice_signal+lambda+1))[y*w+x];
								bg[lambda]+=v;
								bg_sigma[lambda]+=v*v;	
							}
						}
					}
				}
				for (lambda=0; lambda<n; lambda++)
				{
					bg[lambda]/=(float)no_bgpix;
					bg_sigma[lambda]/=(float)no_bgpix;
					bg_sigma[lambda]=(float) Math.sqrt(bg_sigma[lambda]-bg[lambda]*bg[lambda]);
				}
			}
		}
		else if (choice.equals(bgitems[2]))
		{
			goodinput=enterSpectrum(-1);	//enterSpectrum in background mode, i.e. have user enter every channel
		}
		else if (choice.equals(bgitems[3]))		//open a dialog where the user enters a number
		{
			GenericDialog bg_dialog=new GenericDialog("Background selection");	
			bg_dialog.addNumericField("Uniform Background:", 0, 0, 5, "counts");
			bg_dialog.showDialog();
			float bg_uniform=(float) bg_dialog.getNextNumber();	//get a number from the user
			if (bg_uniform>0 && bg_uniform<saturation_threshold)	//sanity check
			{
				goodinput=1;
				for (int i=0; i<n; i++)
				{
					bg[i]=bg_uniform; bg_sigma[i]=1;
				}
			}
			else goodinput=BADVALUE;	//bad input
		}
		if (goodinput==GOODINPUT) return 0;
		else minimumBackground();
		if (goodinput==NOROISELECTED)
		{
			//complain to user!
			IJ.showMessage("No ROI selected! Using minimal signal instead!"); 
			return -1;
		}	
		else if (goodinput==BADVALUE)
		{
			//complain
			IJ.showMessage("Invalid value! Using minimal signal instead!"); 
			return -2;			
		}
		return 0;
	}

	/**
	 * loop over all slices and pixel and determine the minimal signal for each channel
	 */
	private void minimumBackground()
	{
		float min;
		int currentslice_signal;
		for (int lambda=0; lambda<n; lambda++)
		{
			for (int slice=0; slice<tz_dim; slice++)
			{
				currentslice_signal=slice*n+lambda+1;
				min=(float) saturation_threshold;
				for (int pix=0; pix<wh; pix++)
				{	//determine minimum in each channel
					if (((float [])X.getPixels(currentslice_signal))[pix]<min) min=((float [])X.getPixels(currentslice_signal))[pix];
				}
				bg[lambda]=min;
			}
		}
	}

	/**
	 * decide which pixels are saturated or background
	 */
	private void calcBgMask()
	{
		int currentpix, currentslice_signal;
		for (int slice=0; slice<tz_dim; slice++)
		{
			currentpix=slice*wh;
			//initialize the background map with n
			for (int pix=0; pix<wh; pix++)
				bgMask[pix+currentpix]=n;

			//background code: each saturated channel subtracts 2n from bgMask -> a single channel makes it negative
			//each channel below background[lambda]+bg_threshold subtracts 1 from bgMask -> it remains positive unless
			//all channels are below bg
			for (int lambda=0; lambda<n; lambda++)
			{
				currentslice_signal=slice*n+lambda+1;
				for (int pix=0; pix<wh; pix++)
				{	//for each pixel, decide whether it is below background or above saturation
					if (((float [])X.getPixels(currentslice_signal))[pix]<bg[lambda]+no_std*bg_sigma[lambda]+bg_threshold) bgMask[currentpix+pix]--;
					if (((float [])X.getPixels(currentslice_signal))[pix]>saturation_threshold) bgMask[currentpix+pix]-=2*n;
				}
			}
		}
	}


	/**
	 * Prompt user for input, get background and start spectra
	 */
	private int parameterDialog()
	{
		GenericDialog nsources= new GenericDialog("Poisson NMF");
		nsources.addNumericField("Number of Sources", r, 0);
		nsources.showDialog();
		if (nsources.wasCanceled())
		{	
			return CANCELLED;
		}
		else
		{
			r=(int)nsources.getNextNumber();
			if (r<1 || r>10) 
			{
				IJ.error("Bad number of Sources!");
				return -1;
			}
			allocate_spectra();
			//get preferences from IJ preferences file
			getPrefs();
			if (datatype.equals(datatype_choices[1]))
				read_channel_boundaries_LeicaSP2(config_file);
			else if (datatype.equals(datatype_choices[2]))
				read_channel_boundaries_ZeissLSM();
			
			GenericDialog pnmf_dialog=new GenericDialog("Poisson NMF");
			pnmf_dialog.addNumericField("_Number of Iterations", maxit, 0);
			pnmf_dialog.addNumericField("Subsamples", subsamples, 0);
			pnmf_dialog.addNumericField("Segregation Bias", segbias, 0);
			pnmf_dialog.addNumericField("Saturation Threshold", saturation_threshold, 0);
			pnmf_dialog.addNumericField("_Background Threshold", bg_threshold, 0);
			add_spectra_choice(pnmf_dialog);
			pnmf_dialog.addMessage("\n");
			pnmf_dialog.addCheckbox("Specify Spectral Channels?", false);
			pnmf_dialog.setOKLabel("Run!");

			pnmf_dialog.showDialog();
			if (pnmf_dialog.wasCanceled())
			{	
				return CANCELLED;
			}
			else
			{
				//RUNTIME
				//number of iterations
				maxit=(int) pnmf_dialog.getNextNumber();
				if (maxit<1 || maxit>10000) IJ.error("Bad number of Iterations!");
				//Subsampling
				subsamples=(int) pnmf_dialog.getNextNumber();
				if (subsamples<1) subsamples=1;
				//Segregation bias
				segbias=(float) pnmf_dialog.getNextNumber();

				//INITIAL CONDITIONS
				//Saturation Threshold
				saturation_threshold=(float) pnmf_dialog.getNextNumber();
				if (saturation_threshold<0 ) IJ.error("Threshold has to be positive!");
				//Background Threshold
				bg_threshold=(float) pnmf_dialog.getNextNumber();
				if (bg_threshold<0 || bg_threshold>saturation_threshold) IJ.error("Lower threshold has to be negative\n" +
				"and below saturation!");
				determine_channel_order();
				read_spectra_choice(pnmf_dialog);
				//Spectral Windows
				if (pnmf_dialog.getNextBoolean()) 
				{
					calc_channels();
					determine_channel_order();
				}
			}
		}
		setPrefs();
		return 0;
	}


	/**
	 * function that prompts the user for a ROI and calculates the mean spectrum of the 
	 * pixels in that ROI
	 * @param dye
	 */
	private int getSpectrum(int dye)
	{
		boolean goodroi;
		String msg="Select a ROI for dye ";
		msg=msg.concat(Integer.toString(dye+1));
		new WaitForUserDialog("Start Spectra", msg).show();
		Roi spec_roi=imp.getRoi();
		if (spec_roi==null) goodroi=false;
		else goodroi=true;
		if (goodroi)
		{
			//determine time and depth slice currently on display
			int currentslice_signal=imp.getCurrentSlice();
			currentslice_signal--;	//correct for base 1 counting
			currentslice_signal-=currentslice_signal%n;	//determine slice of first channel.
			int npix=0; float v;
			for (int lambda=0; lambda<n; lambda++)
			{
				S[dye][lambda]=0;
			}
			for (int x=0; x<w; x++) //loop over all pixels (terribly inefficient)
			{
				for (int y=0; y<h; y++)
				{
					if (spec_roi.contains(x,y))	//calc mean from pixels in ROI
					{	
						npix++;
						for (int lambda=0; lambda<n; lambda++)
						{
							v=((float [])X.getPixels(currentslice_signal+lambda+1))[y*w+x];
							S[dye][lambda]+=v;
						}
					}
				}
			}
			for (int lambda=0; lambda<n; lambda++)
			{
				S[dye][lambda]/=(float)npix;	//normalize
				S[dye][lambda]-=bg[lambda];		//subtract previouly determined background
			}
			return GOODINPUT;
		}
		else 
		{	//complain since user did not select a ROI
			IJ.showMessage("Select ROI!");
			return NOROISELECTED;
		}
	}

	/**
	 * Open a dialog with a numeric field for each channel such that the user can enter 
	 * the spectrum of the dye 
	 * @param dye
	 * @return
	 */
	private int enterSpectrum(int dye)
	{
		String msg;
		if (dye<0) msg="Enter Background Spectrum";
		else msg="Enter spectrum of dye ";
		msg=msg.concat(Integer.toString(dye+1));
		GenericDialog enter_spectra=new GenericDialog(msg);
		for (int lambda=0; lambda<n; lambda++)
		{	//add numeric field with the current spectrum as default value
			if (dye<0) enter_spectra.addNumericField(Integer.toString(lambda+1), bg[lambda], 4);
			else enter_spectra.addNumericField(Integer.toString(lambda+1), S[dye][lambda], 4);
		}
		enter_spectra.showDialog();
		if (enter_spectra.wasCanceled())
		{	
			return CANCELLED;
		}
		else
		{
			for (int lambda=0; lambda<n; lambda++)
			{	//get user input
				if (dye<0) bg[lambda] =enter_spectra.getNextNumber();
				else	S[dye][lambda]=enter_spectra.getNextNumber();
			}
		}
		return GOODINPUT;
	}

	/**
	 * open dialog that prompts the user for the way the background and the spectra are
	 * to be determined
	 * @return
	 */
	private int add_spectra_choice(GenericDialog choose_spectra)
	{
		//filter files that end on .emn (for emission spectra)
		FilenameFilter spec_file_filter = new FilenameFilter() {
			public boolean accept(File dir, String name) {
				return name.endsWith(".emn");
			}
		};
		//produce a list of files that end on emn in directory /plugins/SpectraLibrary

		//TODO: sort the list of files alphanetically
		//TODO: get the proper path the plugin directory. The current implementation only works if IJ is called from inside its directory
		String curDir = System.getProperty("user.dir");
		curDir=curDir.concat("/plugins/SpectraLibrary");
		File spec_directory=new File(curDir);
		String[] spectra_list=spec_directory.list(spec_file_filter);
		String[] bgitems=new String [4];
		String[] items;
		if (spectra_list!=null)
		{
			items=new String [spectra_list.length + 4];
		}else{
			items=new String [4];			
		}
		//the three choices for the background
		bgitems[0]="Minimal values";
		bgitems[1]="ROI selection";
		bgitems[2]="manually";
		bgitems[3]="flat";
		//the choices for the spectra + the content of the directory /plugins/SpectraLibrary
		items[0]="Gaussian";
		items[1]="ROI selection";
		items[2]="manually";
		items[3]="--------";

		for (int i = 4; i < items.length; i++) {
			items[i]=spectra_list[i-4];
		}

		String label;
		String default_choice=bgitems[0];
		label="_Background spectrum";
		for (int i=0; i<bgitems.length; i++)
		{
			if (bg_choice!=null)
				if (bg_choice.equals(bgitems[i])) default_choice=bgitems[i];
		}
		choose_spectra.addChoice(label, bgitems, default_choice);
		for (int dye = 0; dye < r; dye++) {
			default_choice=items[0];
			for (int i=0; i<items.length; i++)
			{
				if (spec_choice[dye]!=null)
					if (spec_choice[dye].equals(items[i])) 	default_choice=items[i];
			}
			label="Dye_";
			label=label.concat(Integer.toString(dye+1));
			label=label.concat(" initial spectrum");			
			choose_spectra.addChoice(label, items, default_choice);
		}
		String [] labels=new String [r];
		choose_spectra.addMessage("Keep spectra of dyes fixed?");
		for (int dye = 0; dye < r; dye++) 
		{
			labels[dye]="Dye ";
			labels[dye]=labels[dye].concat(Integer.toString(dye+1));
		}
		choose_spectra.addCheckboxGroup(1, r, labels, spectra_fixed);
		return 0;

	}

	/**
	 * handle the choice of spectra from the parameter dialog
	 *
	 */
	private int read_spectra_choice(GenericDialog choose_spectra)
	{
		//the three choices for the background
		String[] items=new String [4];
		String curDir = System.getProperty("user.dir");
		curDir=curDir.concat("/plugins/SpectraLibrary");

		//the choices for the spectra + the content of the directory /plugins/SpectraLibrary
		items[0]="Gaussian";
		items[1]="ROI selection";
		items[2]="manually";
		items[3]="--------";
		//deal with background first, it is needed for the determination of the spectra later
		bg_choice=choose_spectra.getNextChoice();
		getBackground(bg_choice);
		calcBgMask();

		gauss_spectra();						//initialize with gaussians
		//get choice for each spectrum and assign spectra
		for (int dye = 0; dye < r; dye++) {
			spec_choice[dye]=choose_spectra.getNextChoice();
			if (spec_choice[dye].equals(items[0]) || spec_choice[dye].equals(items[3])) {}	//Gaussian
			else if (spec_choice[dye].equals(items[1])) getSpectrum(dye);			//ROI
			else if (spec_choice[dye].equals(items[2])) enterSpectrum(dye);		//Manually
			else readSpectrum(dye, curDir, spec_choice[dye]);				//from file
		}
		//read whether spectra are to be kept fixed
		for (int dye = 0; dye < r; dye++) {
			spectra_fixed[dye]=choose_spectra.getNextBoolean();
		}
		
		//Normalize the spectra
		for (int z=0;z<r;z++) {				//for each dye, calculate the 'E'-norm of the spectrum
			power[z]=0.0;
			for (int y=0;y<n;y++)
				power[z]+=Math.pow(S[z][y],E);
			power[z]=(float) Math.pow(power[z], 1f/E);
		}
		for (int z=0;z<r;z++) {			//rescale S
			for (int y=0;y<n;y++)
			{
				S[z][y]/=power[z];
			}
		}
		//save initial spectra for nostalgia and normalize them for display
		int ii;
		for (int dye = 0; dye < r; dye++) {
			for (int lambda=0; lambda<n; lambda++)	
				{
					ii=channel_order[lambda];
					initialS[dye][lambda]=(float) S[dye][ii];
				}
		}
		float maxS=0;
		for (int z=0;z<r;z++) {
			for (int y=0;y<n;y++) {
				ii=channel_order[y];
				if (initialS[z][y]/(channel_lambdas[ii][1]-channel_lambdas[ii][0])>maxS) 
					maxS=(float) ((float) initialS[z][y]/(channel_lambdas[ii][1]-channel_lambdas[ii][0]));
			}
		}
		for (int z=0;z<r;z++) {
			for (int y=0;y<n;y++) {
				ii=channel_order[y];
				initialS[z][y]=(float) (initialS[z][y]/(channel_lambdas[ii][1]-channel_lambdas[ii][0])/maxS);
			}
		}
		return 0;
	}

	/**
	 * Read spectrum form file. 
	 * @param dye
	 * @param dir
	 * @param spec
	 */
	private void readSpectrum(int dye, String dir, String spec)
	{
		try {
			BufferedReader spectrum=new BufferedReader(new FileReader(dir+"/"+spec));
			String line;
			StringTokenizer entries;
			Vector<Double> lambda=new Vector<Double>();
			Vector<Double> emission=new Vector<Double>();
			//read in all line from the file and assume the first column is wavelength
			//the second column is spectrum
			while((line=spectrum.readLine())!=null)
			{
				entries=new StringTokenizer(line);
				if (entries.countTokens()>1){
					lambda.addElement(Double.valueOf(entries.nextToken()));
					emission.addElement(Double.valueOf(entries.nextToken()));
				}
			}
			spectrum.close();
			//Integrate the spectrum read above for each spectral channel
			int c=0;
			double step, result;	//integration step
			for (c = 0; c < n; c++) 
			{
				step=(channel_lambdas[c][1]-channel_lambdas[c][0])/20;	//integration step
				result=0;
				for (double lambdax=channel_lambdas[c][0];lambdax<channel_lambdas[c][1]; lambdax+=step)
				{	//use midpoint of interval as pivot for interpolated value of spectrum
					result+=interpolate(lambda, emission, lambdax+0.5*step)*step;
				}
				S[dye][c]=result;
			}

		} catch (FileNotFoundException e) {
			IJ.error("File not found");
		} catch (IOException e) {
			IJ.error("Bad file!");
		}

	}

	/**
	 * interpolate the list of lambda and spectral values read from the file linearly. 
	 * @param lambda
	 * @param F
	 * @param lambda_0
	 * @return
	 */
	private double interpolate(Vector<Double> lambda, Vector<Double> F, double lambda_0)
	{
		int i=0;
		double result=0;
		if (lambda_0<lambda.get(0)) //value outside the interval -> extrapolate
		{
			double lambda1=lambda.get(0), lambda2=lambda.get(1), F1=F.get(0), F2=F.get(1);
			result=F1+(F2-F1)/(lambda2-lambda1)*(lambda_0-lambda1);
			result=Math.min(F1, F2);
		}
		else if (lambda_0>lambda.get(lambda.size()-1))  //value outside the interval -> extrapolate
		{
			double lambda1=lambda.get(lambda.size()-2), lambda2=lambda.get(lambda.size()-1), F1=F.get(lambda.size()-2), F2=F.get(lambda.size()-1);
			result=F1+(F2-F1)/(lambda2-lambda1)*(lambda_0-lambda1);
			//result=F2;
		}
		else //interpolate linearly
		{
			//find the relevant interval
			while (i<lambda.size() && lambda.get(i)<lambda_0)	i++;
			if (F.size()>i){
				double lambda1=lambda.get(i-1), lambda2=lambda.get(i), F1=F.get(i-1), F2=F.get(i);
				result=F1+(F2-F1)/(lambda2-lambda1)*(lambda_0-lambda1);
				//result=(lambda_0-lambda1<lambda2-lambda_0)?F1:F2;
			}else result=0;
		}
		return Math.max(result,0);	//make sure the value is positive
	}

	/**
	 * calculate the center wavelength of all channels -- assumes the minimal and maximal
	 * wavelength entered by the user are center wavelength of first and last channel
	 * @param lambda_min
	 * @param lambda_max
	 * @return
	 */
	private int calc_channels()
	{
		String msg="Enter emission channels.";
		GenericDialog enter_channels=new GenericDialog(msg);
		for (int lambda=0; lambda<n; lambda++)
		{	//add numeric field with the current spectrum as default value
			enter_channels.addNumericField("Channel "+Integer.toString(lambda+1)+" lower boundary: ", channel_lambdas[lambda][0], 6);
			enter_channels.addNumericField("Channel "+Integer.toString(lambda+1)+" upper boundary: ", channel_lambdas[lambda][1], 6);
		}
		enter_channels.showDialog();
		if (enter_channels.wasCanceled())
		{	
			return -1;
		}
		else
		{
			for (int lambda=0; lambda<n; lambda++)
			{	//get user input
				channel_lambdas[lambda][0]=enter_channels.getNextNumber();
				channel_lambdas[lambda][1]=enter_channels.getNextNumber();
			}
		}
		return 0;
	}

	/**
	 * Do the update for the spectra
	 */
	private void updateS()
	{
		int x; int y; int z;
		float[] sumRowA = new float[r];
		float f;

		// calc estimated signal A*S
		for (x=0;x<noPix;x++)
			for (y=0;y<n;y++)					
				for (z=0,ASsub[x][y]=0;z<r;z++)
					ASsub[x][y]+=Asub[x][z]*S[z][y];
		//calculate normalizer for the update rule
		for (z=0;z<r;z++)
			for (x=0,sumRowA[z]=0;x<noPix;x++)
				sumRowA[z]+=Asub[x][z];
		//do update
		for (z=0;z<r;z++) {
			if (spectra_fixed[z]==false)
			{
				for (y=0;y<n;y++) {
					f=0;
					for (x=0;x<noPix;x++) 
						f+=Asub[x][z]*Xsub[x][y]/ASsub[x][y];
					S[z][y]*=f/sumRowA[z];				
				}
			}
		}
	}

	/**
	 * Do the update for the concentrations
	 */
	private void updateA()
	{
		int x; int y; int z;
		double f;
		float[] sumColS = new float[r];
		double onenorm, twonorm;

		// calc estimated signal A*S
		for (x=0;x<noPix;x++)
			for (y=0;y<n;y++)					
				for (z=0,ASsub[x][y]=0;z<r;z++)
					ASsub[x][y]+=Asub[x][z]*S[z][y];

		if (segbias>0)		//calculate the segregation bias term only if needed
		{
			for (x=0;x<noPix;x++)
			{
				for (z=0, onenorm=0, twonorm=0;z<r;z++)
				{
					f=Asub[x][z];
					onenorm+=f;				// sum x_i
					twonorm+=f*f;			// sum x_i^2
				}
				onenorm/=twonorm;			
				twonorm=Math.sqrt(twonorm);	//twonorm is now the proper two norm sqrt(sum(x_i^2))
				onenorm/=twonorm;			//onenorm is now sum(x_i)/(sum x_i^2)^(3/2)
				for (z=0;z<r;z++)
				{							//calculate the gradient of |x|/||x|| and multiply by segbias
					Segbiasterm[x][z]=segbias*(onenorm*Asub[x][z]-1.0/twonorm);
				}		
			}
		}

		//calculate normalizer for update rule
		for (z=0;z<r;z++)
			for (y=0,sumColS[z]=0;y<n;y++)
				sumColS[z]+=S[z][y];
		//do update
		for (z=0;z<r;z++) {
			for (x=0;x<noPix;x++) {
				f=0;
				for (y=0;y<n;y++) 
					f+=Xsub[x][y]/ASsub[x][y]*S[z][y];
				Asub[x][z] *= (f+Segbiasterm[x][z])/sumColS[z];
				if (Asub[x][z]<conc_nothing) Asub[x][z]=conc_nothing;	//make sure nothing becomes negative
			}
		}
	}

	/**
	 * function that calculates the mid_wavelength of each channel and determines 
	 * the order of the channels with increasing wavelength
	 * this is important to produce reasonable spectra plots if the images in the source stack are
	 * not ordered, as it can happen with Leica SP2 data.
	 */
	private void determine_channel_order()
	{
		double [] mid_wavelength=new double [n];
		double [] mid_wavelength_sorted=new double [n];
		for (int channel=0; channel<n; channel++) {
			mid_wavelength[channel]=(channel_lambdas[channel][0]+channel_lambdas[channel][1])/2;
			mid_wavelength_sorted[channel]=(channel_lambdas[channel][0]+channel_lambdas[channel][1])/2;			
		}
		Arrays.sort(mid_wavelength);
//		double prev_min=-1, min;
//		int c=0;
//		for (int channel=0; channel<n; channel++) 
//		{
//			min=1000;
//			for(int channel2=0; channel2<n; channel2++)
//			{
//				if (mid_wavelength[channel2]<min && mid_wavelength[channel2]>prev_min)
//				{min=mid_wavelength[channel2]; c=channel2;}
//			}
//			channel_order[channel]=c;
//			prev_min=min;
//		}
		for (int channel=0; channel<n; channel++) 
		{
			for(int channel2=0; channel2<n; channel2++)
			{
				if (mid_wavelength[channel2]==mid_wavelength_sorted[channel])
					channel_order[channel2]=channel;
			}
		}
		for (int channel=0; channel<n; channel++) 
		{
			inverse_channel_order[channel_order[channel]]=channel;
		}
	}
	
	/**
	 * Function that calculates the pseudo-inverse of the spectra matrix
	 */
	private void calc_pinv(boolean mod)
	{
		pinvS=new double [r][n];
		double [][] STS=new double [r][r];
		double [][] STSch=new double [r][r];
		double [] y=new double [r];

		for (int i = 0; i < STS.length; i++) {
			for (int j = 0; j < STS[i].length; j++) {
				STS[i][j]=0;
				for (int k = 0; k < n; k++) {
					if (mod)
						STS[i][j]+=modifiedS[i][k]*modifiedS[j][k];
					else
						STS[i][j]+=S[i][k]*S[j][k];						
				}
			}
		}

		double sum,diff;
		int i,j,l, k1, k2;
		//cholesky decomposition
		for (i = 0; i < r; i++) {
			sum = 0.0;
			for (j = 0; j < i; j++) {
				sum += STSch[i][j]*STSch[i][j];
			}
			diff = STS[i][i] - sum;
			STSch[i][i] = Math.sqrt(diff);
			for (l = i+1; l < r; l++) {
				sum = 0.0;
				for (j = 0; j < i; j++) {
					sum += STSch[l][j]*STSch[i][j];
				}
				STSch[l][i] = (STS[l][i] - sum)/STSch[i][i];
				STSch[i][l]=0;
			}
		}
		//solving the lower diagonal system
		for (k1 = 0; k1 < r; k1++) {
			for (k2 = 0; k2 < r; k2++) y[k2]=0;
			y[k1]=1;
			for (i = 0; i < r; i++) {
				sum = 0.0;
				for (j = 0; j < i; j++) {
					sum += STSch[i][j]*STS[j][k1];
				}
				STS[i][k1] = (y[i] - sum)/STSch[i][i];
			}	
		}
		//calculate the product of the inverse of the lower diagonal
		for (k1 = 0; k1 < r; k1++) {
			for (k2 = k1; k2 <r; k2++) {
				STSch[k1][k2] = 0.0;
				for (j = k2; j <r; j++) {
					STSch[k1][k2] += STS[j][k1]*STS[j][k2];
				}
				STSch[k2][k1]=STSch[k1][k2];
			}	
		}
		//calculate the pseudo inverse
		for (k1 = 0; k1 < r; k1++) {
			for (k2 = 0; k2 < n; k2++) {
				pinvS[k1][k2] = 0.0;
				for (j = 0; j < r; j++) {
					pinvS[k1][k2] += STSch[k1][j]*S[j][k2];
				}
			}	
		}
	}
	
	private void results_panel()
	{
		PoissonNMFresultsPanel results=new PoissonNMFresultsPanel(this);
	}

	/**
	 * Open a series of Leica Images and merge them into a stack
	 * @return
	 */
	private int openLeicaSP2()
	{
		//filter files that end on .tif
		FilenameFilter im_file_filter = new FilenameFilter() {
			public boolean accept(File dir, String name) {
				return name.endsWith(".tif");
			}
		};
		//filter files that end on .txt
		FilenameFilter txt_file_filter = new FilenameFilter() {
			public boolean accept(File dir, String name) {
				return name.endsWith(".txt");
			}
		};
		
		OpenDialog od=new OpenDialog("Select one image of the series", null);
		String LeicaFile=od.getFileName();
		//Check sanity of input or for cancellation
		if (LeicaFile==null) return CANCELLED;
		else if (LeicaFile.endsWith(".txt") || LeicaFile.endsWith(".lei"))
		{
			IJ.showMessage("Select one image of the desired series, not the *.txt or *.lei file");
			return openLeicaSP2();
		}
		String series;
		try{
			series=LeicaFile.substring(LeicaFile.indexOf("Series"),LeicaFile.indexOf("Series")+9);
		}catch (Exception e)
		{
			IJ.showMessage("Unable to extract series from file name");
			return CANCELLED;
		}
		// read all images that of the same series and add them to a newly generated stack
		String LeicaDir = od.getDirectory();
		File dir=new File(LeicaDir);
		String[] im_list=dir.list(im_file_filter);
		ImageStack datastack=null;
	    for (int i=0; i<im_list.length; i++) {
	        IJ.showProgress(i, im_list.length);
	        if (im_list[i].contains(series)) {
	            ImagePlus im=IJ.openImage(LeicaDir+im_list[i]);
	            if (imp==null) {
	                w=im.getWidth(); h=im.getHeight();
	                imp=new ImagePlus();
	                datastack=new ImageStack(w,h);
	                datastack.addSlice("Channel "+i, im.getProcessor());
	            } else {
	                datastack.addSlice("Channel "+i, im.getProcessor());
	            }
	        }
	    }
	    //Display the stack and remember the data configuration file, assuming it is the only
	    //text file in the directory
		imp = new ImagePlus("Data", datastack);
		imp.show();
		String[] txt_list=dir.list(txt_file_filter);
		config_file=LeicaDir+txt_list[0];
		return 0;
	}

	/**
	 * open a regular image stack using the imageJ command. it also opens the Zeiss LSM
	 * stacks.
	 * @return
	 */
	private int openStack()
	{
		IJ.run("Open...");
		imp=WindowManager.getCurrentImage();
		return 0;
	}

	/**
	 * Function that parses the Leica SP2 configuration file
	 * @param conf_file
	 */
	private void read_channel_boundaries_LeicaSP2(String conf_file)
	{
		try {
			BufferedReader LeicaInfoFile=new BufferedReader(new FileReader(conf_file));
			String line;
			StringTokenizer entries;
			//read in all line from the file and assume the first column is wavelength
			//the second column is spectrum
			int channel=0;
			int sequence=1;
			String seq_base="HARDWARE PARAMETER #0 SEQUENCE ";
			String scan_base="SCANNER INFORMATION #0";
			//read everything in the file until the parameters for the first data sequence are found 
			while((line=LeicaInfoFile.readLine())!=null && (!line.contains(seq_base+sequence) && !line.contains(scan_base)))
			{}
			if (line==null) return;
			else if (line.contains(seq_base+sequence))
			{	
				while((line=LeicaInfoFile.readLine())!=null)
				{
					if (line.contains("SP Mirror"))
					{
						entries=new StringTokenizer(line);
						entries.nextToken();
						entries.nextToken();
						entries.nextToken();
						String lr=entries.nextToken();
						String wavelength=entries.nextToken();
						double wl=0; 
						//read entry where a number is expected
						try{wl=Double.valueOf(wavelength);}
						catch (Exception e) {;}
						//save it as a channel boundary according to label left and right
						if (lr.equals("(left)")) channel_lambdas[channel][0]=wl;
						else if (lr.equals("(right)")) {channel_lambdas[channel][1]=wl; channel++;}
						//when reading the right boundary, increment the channel index
						//jump out of function once "enough" channel boundaries have been read.
						if (channel>=n) return;
					}
					if (line==seq_base+sequence+1) sequence++;
				}
			}
			else if (line.contains(scan_base)){
				while((line=LeicaInfoFile.readLine())!=null)
				{
					if (line.contains("[Lambda"))
					{
						entries=new StringTokenizer(line);
						String wavelength=entries.nextToken();
						wavelength=entries.nextToken();
						wavelength=entries.nextToken();
						double wl=0; 
						//read entry where a number is expected
						try{wl=Double.valueOf(wavelength);}
						catch (Exception e) {;}
						if (line.contains("[LambdaBeginLeft")) channel_lambdas[0][0]=wl;
						else if (line.contains("[LambdaBeginRight")) channel_lambdas[0][1]=wl;
						else if (line.contains("[LambdaEndLeft")) channel_lambdas[n-1][0]=wl;
						else if (line.contains("[LambdaEndRight")) channel_lambdas[n-1][1]=wl;
					}
				}
				for (channel=1; channel<n-1; channel++)
				{
					channel_lambdas[channel][0]=channel_lambdas[0][0]+(channel_lambdas[n-1][0]-channel_lambdas[0][0])/(n-1)*channel;
					channel_lambdas[channel][1]=channel_lambdas[0][1]+(channel_lambdas[n-1][1]-channel_lambdas[0][1])/(n-1)*channel;
				}				
			}
		} catch (FileNotFoundException e) {
			return;
		} catch (IOException e) {
			return;
		}
	}

	/**
	 * mid wavelength of the channels are included in the slice labels of the image stacks
	 * opened from LSM files.
	 */
	private void read_channel_boundaries_ZeissLSM()
	{
		String [] slicelabels=X.getSliceLabels();
		double [] mid_wavelength=new double [n];
		boolean goodboundaries=true;
		if ((imp.getOriginalFileInfo().fileName).endsWith(".lsm"))
		{
			for (int i = 0; i < mid_wavelength.length; i++) {
				mid_wavelength[i]=Double.valueOf(slicelabels[i]);
			}
			//if the result is reasonable, use the information to generate wavelength filling
			//emission windows
			if (mid_wavelength[0]>100 && mid_wavelength[0]<2000)
				for (int i = 0; i < mid_wavelength.length-1; i++) {
				{
					channel_width=mid_wavelength[i+1]-mid_wavelength[i];
					channel_lambdas[i][0]=mid_wavelength[i]-0.5*channel_width;
					channel_lambdas[i][1]=mid_wavelength[i]+0.5*channel_width;
				}
				channel_lambdas[n-1][0]=mid_wavelength[n-1]-0.5*channel_width;
				channel_lambdas[n-1][1]=mid_wavelength[n-1]+0.5*channel_width;
			}else {
				goodboundaries=false;
				IJ.showMessage("Channel boundaries don't make sense, using equal spacing instead!");
			}
		}
		else
		{
			IJ.showMessage("Expected Zeiss LSM file. Using equally space channel boundaries instead!");
			goodboundaries=false;
		}
		// If we had trouble up there, use some ad hoc spacing.
		if (goodboundaries==false)
		{
			double start=500;
			channel_width=11;
			for (int i = 0; i < n; i++) {
				channel_lambdas[i][0]=start+i*channel_width;
				channel_lambdas[i][1]=start+(i+1)*channel_width;
			}
		}
	}

}

