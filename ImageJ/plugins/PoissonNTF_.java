import ij.*;
import ij.process.*;
import ij.gui.*;
import ij.io.*;

import java.awt.*;
import java.io.*;
import java.util.StringTokenizer;
import java.util.Vector;

import ij.plugin.filter.*;

/**
 * PoissonNTF_ 
 * imagej plugin for poissonNTF
 *
 * @author <a href="mailto:fabian.theis@helmholtz-muenchen.de">Fabian Theis</a>
 * @author <a href="mailto:neher@kitp.ucsb.edu">Richard Neher</a>
 */
public class PoissonNTF_ implements PlugInFilter {
	String PoissonNTFversion="0.3.1";

	ImagePlus [] imp=null;	//input, holds the stack X
	ImageStack [] X=null; 	//emission signal
	ImagePlus img=null;	//output, hold the stack A
	ImageStack A=null;	//concentrations 

	//dimension of the images
	int ndyes=3;
	int nemn;
	int nexc;
	int w; //width;
	int h; //height;
	int wh; //product
	int tz_dim; //number of product of dimensions 4 and 5
	int totalchannels;
	int [][] dim; //Dimensions of each data stack;
	int [] nchannels;
	int [][] channeltoemission;
	int []lowerexc=null;
	int []upperexc=null;
	
	//internal error codes
	int CANCELLED=-1;
	int NOROISELECTED=-2;
	int GOODINPUT=1;
	int BADVALUE=-3;

	// Strings that contain the current normalization state
	String concNorm = "";
	String emnNorm = "";	
	String excNorm = "";
	int currentnormalization=0;
	
	//matrices that hold the NTF objects and intermediates
	int []bgMask;
	double[][][] Xsub;	//signal
	double[][] Asub;	//concentration
	double[][][] QASsub;	//expected concentration
	double[][] S;		//Emission spectra
	double[][] Q;		//Excitation spectra
	float[][] initialS;	//initial emission spectra
	double[][] pinvS;	//pseudo inverse of unfolded set of equations

	//auxillary quantities and parameters
	Vector <Vector<Double>> ChannelsL=null;
	Vector <Vector<Double>> ChannelsR=null;
	Vector <String> excitation_wavelength=null;
	Vector <Double> emission_wavelength=null;
	Vector <Double> channel_width=null;
	int [] channel_order;
	int [] inverse_channel_order;
	double[][] Segbiasterm;		//temporary variable to store gradient contribution of the seg bias
	double[][] bg; 		//background signal (vector, one entry for each spectral channel)
	double[][] bg_sigma;	//std of the above
	double [][][] Sexc;		//array holding the spectra from each excitation separately [nexc][nemm][ndyes]
	double[] power; 	//buffer
	float signal_nothing=1;		//positive but irrelevant signal strength
	double conc_nothing=0.001;		//positive but irrelevant concentration
	double spec_nothing=0.001;		//positive but irrelevant emission
	double E=1; 		//E=1 one-norm for spectra, E=2 two norm
	String excNormType="Max";	//prescription to normalize the excitation spectra
	double segbias=0;	//weight of the segregation bias
	int noPix;			//number of pixels in the current subsample
	boolean isHyperstack;
	double bg_threshold=20;		// Threshold, below which data is not used  
	double no_std=2;				// number of standard deviations of background noise, below which data is not used
	double saturation_threshold=4000;	//threshold above which the data is discarded
	boolean []spectra_fixed;			//vector of booleans to decide which spectra to fix during optimization
	boolean [] imageselected;
	int maxit=100;						//number of iterations
	int subsamples=3;					//subsample levels, each with 10fold less data
	int nmf_preiterations=50;				// number of iterations on concentrations without spectra update
	int no_preiterations=20;				// number of iterations on concentrations without spectra update
	int no_postiterations=5;				// number of iterations on concentrations without spectra update before showing concentrations
	boolean cancelled=false;
	boolean interrupted=false;
	String [] spec_choice;			//array of strings with the choice of initial spectra
	String bg_choice;				//choice for the background mode
	PlotWindow plotEmission=null;			//window handle to display the spectra
	PlotWindow plotExcitation=null;			//window handle to display the spectra
	String datatype;
	String [] datatype_choices;
	String config_file;
	/**
	 * initalization of the routine, basically only loading of data
	 * @param arg
	 * @param imp
	 * @return
	 * @see ij.plugin.filter.PlugInFilter#setup(java.lang.String, ij.ImagePlus)
	 */
	public int setup(String arg, ImagePlus impin) {
		datatype_choices=new String [3];
		datatype_choices[0]="Regular stacks";
		datatype_choices[1]="Leica SP2";
		datatype_choices[2]="Zeiss LSM";
		datatype=datatype_choices[1];
		if (imp==null)	//if no image stack is passed, open dialog that lets the user choose the data
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
				if (datatype.equals(datatype_choices[1]))	//call routine to open Leica Data
				{if (openLeicaSP2()==CANCELLED) return DONE;}
				else if (datatype.equals(datatype_choices[2]))	//any other choice, prompt the user for stacks (Zeiss is not yet implemented)
				{if (openStacks()==CANCELLED) return DONE;}
				else if (datatype.equals(datatype_choices[0]))	////any other choice, prompt the user for stacks
				{if (openStacks()==CANCELLED) return DONE;}
				else
				{;}
			}
		}
		if (imp!=null)	//if a stack was open, use this data
		{
			//get image dimensions and check whether there are several image slices
			for (int exc=0; exc<nexc; exc++)
				dim[exc]=imp[exc].getDimensions();


			int d=1; for (int i = 2; i < dim[0].length; i++) d*=dim[0][i];
			if (d==1)
			{
				IJ.error("Stack required!");
				return DONE; 	
			}
			//if stack is ok, convert to 32 gray scale and assign dimensions.
			isHyperstack=imp[0].isHyperStack();
			for (int exc=0; exc<nexc; exc++)
			{
				if ((imp[exc].getBytesPerPixel()<32)) new StackConverter(imp[exc]).convertToGray32();
				X[exc] = imp[exc].getStack();
			}
			w = dim[0][0]; //X.getWidth();
			h = dim[0][1]; //X.getHeight();
			wh = w*h;	//number of pixels per frame
			//for hyperstacks, save number of frames in tz_dim
			if (isHyperstack) {tz_dim=dim[0][3]*dim[0][4];}
			else {tz_dim=1;}
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
		//function exits and proceeds with run()
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
			IJ.error("PoissonNTF requires ImageJ version 1.39t or greater!");
			return;
		}
		//user input of basic parameters - opens dialog
		if (parameterDialog()<0) return;
		int []iterations=new int [subsamples];
		double sum=0;
		/**
		 * normalize the segregation bias such that the contribution of the additional term to the cost
		 * function is of the order of that of a background pixel when segbias was 1
		 * the necessary magnitude of the segregation bias  depends on the overlap of the 
		 * spectra and the distribution of dye concentrations
		 */
		segbias*=nexc*0.5*Math.log(Math.max(bg_threshold,100));
		//calculate the number of iterations made for each subsample level
		for (int s=subsamples-1; s>=0; s--) sum+=Math.pow(2,s);
		for (int s=subsamples-1; s>=0; s--) iterations[s]=(int) (maxit*Math.pow(2,s)/sum);
		double sum2=0;
		for (int s=subsamples-1; s>=0; s--) sum2+=Math.pow(0.1, s)*iterations[s];

		//the iteration loop -- the number of iterations is specified by the user
		double count=0;
		PoissonNTFprogress progress=new PoissonNTFprogress(this);

		subsample(Math.pow(0.1, subsamples-1));	//produce a new subsample of size 10^-s
		for (int s=subsamples-1; s>=0; s--)
		{
			subsample(Math.pow(0.1, s));	//produce a new subsample of size 10^-s
			if (s==subsamples-1)			//initialize only for first subsample, otherwise solve for concentration
			{		
				initialization_QAS();
			}
			else solveforA();				//initialize the concentrations with least squares solution
			for (int it=0; it<iterations[s]; it++) {
				IJ.showProgress(count/sum2);
				normalizeNTF();		//normalize spectra and concentrations
				updateS_NTF();		//update the spectra
				updateA_NTF();		//update the concentrations
				updateQ_NTF();		//update the concentrations
				count+=Math.pow(0.1, s);	
				if (cancelled)	return;	//check state of cancel button
				else if (interrupted) { s=-1;  break;}	//break out of loops if interrupted
				plotEmissionSpectra(false);
				plotExcitationSpectra(false);
			}
		}
		if (!(cancelled || interrupted)) progress.closeWindow();
		IJ.showProgress(1);

		//display results
		showConcentrations(true);
		plotExcitationSpectra(true);
		plotEmissionSpectra(true);
		PoissonNTFresultsPanel resultspanel=new PoissonNTFresultsPanel(this);
	}

	/**
	 * To initialize, start with random concentrations but guessed (or specified) spectra
	 * and adjust the concentrations doing NMF
	 */
	private void initialization_QAS() {
		int exc;
		initA();
		for (exc=0; exc<nexc; exc++){	//at first, do iterations on concentrations only, since spectra might be set
			for (int it=0; it<no_preiterations; it++){
				normalizeNMF(exc);
				updateA_NMF(exc);
			}
			for (int it=0; it<nmf_preiterations; it++){
				normalizeNMF(exc);	//do full NMF iterations
				updateS_NMF(exc);
				updateA_NMF(exc);
			}
		}
		splice_spectra();		//since not all illuminations cover all channels, spectra have to be combined.
								//at the same time, the approximate excitations efficiencies are determined
	}

	/**
	 * init Asub with random number
	 */
	private void initA() {
		int x;int dye;

		for (x=0;x<noPix;x++)		//loop over pixels and assign concentrations with random numbers
			for (dye=0;dye<ndyes;dye++) {		//loop over dyes
				Asub[x][dye]= (float) (Math.random()/2+.5);
			}
	}


	/**
	 * choose subsample from above background data points, subsample fraction is frac
	 */
	private void subsample(double frac)
	{
		int exc, emn;
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
				if (bgMask[pix]>0 && Math.random()<frac) {bgMask[pix]=2*totalchannels;}
				if (bgMask[pix]==2*totalchannels) noPix++;
			}
		}
		Xsub=new double [nexc][noPix][nemn];
		Asub=new double[noPix][ndyes];
		Segbiasterm=new double [noPix][ndyes];
		QASsub=new double [nexc][noPix][nemn];
		//loop over slices and pixels and copy the signal into from the ImagePlus
		//into Xsub, pixels are chosen on the basis of bgMask, count is the continuous index
		for (int slice=0; slice<tz_dim; slice++)
		{
			currentpix=slice*wh;
			for (int pix=0; pix<wh; pix++)
			{
				if (bgMask[currentpix+pix]==2*totalchannels)	//if pixel is selected, add it to the data set
				{
					for(exc=0; exc<nexc; exc++){
						currentslice_signal=slice*nchannels[exc];
						for (emn=0; emn<nchannels[exc]; emn++)
						{	//subtract the background signal, set negative pixels to signal_nothing
							Xsub[exc][count][channeltoemission[exc][emn]]= ((float [])X[exc].getPixels(currentslice_signal+emn+1))[pix]-bg[exc][emn];
							if (Xsub[exc][count][channeltoemission[exc][emn]]<signal_nothing) Xsub[exc][count][channeltoemission[exc][emn]]=signal_nothing;
						}
					}
					count++;
				}
			}
		}
	}

	/**
	 * procedure that calcuates the signal using the current model
	 * i.e. the product out of Q,A,S
	 * This is needed in the update routines for NTF
	 */
	private void calc_QASsub()
	{
		int exc, emn, x, dye;
		for(exc=0; exc<nexc; exc++){
			for (x=0;x<noPix;x++){
				for (emn=0;emn<nemn;emn++){
					QASsub[exc][x][emn]=0;
					for (dye=0;dye<ndyes;dye++){
						QASsub[exc][x][emn]+=Q[exc][dye]*Asub[x][dye]*S[emn][dye];
					}}}}
	}

	/**
	 * Do the update for the spectra
	 */
	private void updateS_NMF(int exc)
	{
		int x; int emn; int dye;
		float[] sumRowA = new float[ndyes];
		float f;
		// calc estimated signal Q*A*S for the subspectra
		for (x=0;x<noPix;x++){
			for (emn=0;emn<nchannels[exc];emn++){
				QASsub[exc][x][emn]=0;
				for (dye=0;dye<ndyes;dye++){
					QASsub[exc][x][emn]+=Q[exc][dye]*Asub[x][dye]*Sexc[exc][emn][dye];
				}}}

		//calculate normalizer for the update rule
		for (dye=0;dye<ndyes;dye++){
			for (x=0,sumRowA[dye]=0;x<noPix;x++){
				sumRowA[dye]+=Q[exc][dye]*Asub[x][dye];
			}}
		//do update
		for (dye=0;dye<ndyes;dye++) {
			if (spectra_fixed[dye]==false)
			{
				for (emn=0;emn<nchannels[exc];emn++) {
					f=0;
					for (x=0;x<noPix;x++) {
						f+=Asub[x][dye]*Xsub[exc][x][channeltoemission[exc][emn]]/QASsub[exc][x][emn];
					}
					Sexc[exc][emn][dye]*=f/sumRowA[dye];				
				}
			}
		}
	}

	/**
	 * Do the update for the concentrations in the NMF scheme
	 */
	private void updateA_NMF(int exc)
	{
		int x; int emn; int dye;
		double f;
		float[] sumColS = new float[ndyes];
		double onenorm, twonorm;

		// calc estimated signal A*S
		for (x=0;x<noPix;x++){
			for (emn=0;emn<nchannels[exc];emn++){
				QASsub[exc][x][emn]=0;
				for (dye=0;dye<ndyes;dye++){
					QASsub[exc][x][emn]+=Q[exc][dye]*Asub[x][dye]*Sexc[exc][emn][dye];
				}}}

		if (segbias>0)		//calculate the segregation bias term only if needed
		{
			for (x=0;x<noPix;x++)
			{
				for (dye=0, onenorm=0, twonorm=0;dye<ndyes;dye++)
				{
					f=Asub[x][dye];
					onenorm+=f;				// sum x_i
					twonorm+=f*f;			// sum x_i^2
				}
				onenorm/=twonorm;			
				twonorm=Math.sqrt(twonorm);	//twonorm is now the proper two norm sqrt(sum(x_i^2))
				onenorm/=twonorm;			//onenorm is now sum(x_i)/(sum x_i^2)^(3/2)
				for (dye=0;dye<ndyes;dye++)
				{							//calculate the gradient of |x|/||x|| and multiply by segbias
					Segbiasterm[x][dye]=segbias*(onenorm*Asub[x][dye]-1.0/twonorm);
				}		
			}
		}

		//calculate normalizer for update rule
		for (dye=0;dye<ndyes;dye++){
			sumColS[dye]=0;
			for (emn=0;emn<nchannels[exc];emn++){
				sumColS[dye]+=Q[exc][dye]*Sexc[exc][emn][dye];
			}}
		//do update
		for (dye=0;dye<ndyes;dye++) {
			for (x=0;x<noPix;x++) {
				f=0;
				for (emn=0;emn<nchannels[exc];emn++) {
					f+=Xsub[exc][x][channeltoemission[exc][emn]]/QASsub[exc][x][emn]*Sexc[exc][emn][dye];
				}
				Asub[x][dye] *= (f+Segbiasterm[x][dye])/sumColS[dye];
				if (Asub[x][dye]<conc_nothing) Asub[x][dye]=conc_nothing;	//make sure nothing becomes negative
			}
		}
	}


	/**
	 * Do the update for the emission spectra
	 */
	private void updateS_NTF()
	{
		int x; int emn; int dye, exc;
		float[][] sumRowQA = new float[nemn][ndyes];
		float f;
		// calc estimated signal Q*A*S for the subspectra
		calc_QASsub();

		//calculate normalizer for the update rule
		//note that this normalizer is degenerate. it does not depend on emn
		//the only reason to include it is the implicit dependence on emn through the excitations
		//a particular channel is illuminated
		for (dye=0;dye<ndyes;dye++){
			for (emn=0;emn<nemn;emn++) {
				sumRowQA[emn][dye]=0;
				for(exc=lowerexc[emn]; exc<upperexc[emn]; exc++){
					for (x=0;x<noPix;x++){
						sumRowQA[emn][dye]+=Q[exc][dye]*Asub[x][dye];
					}}}}
		//do update
		for (dye=0;dye<ndyes;dye++) {
			if (spectra_fixed[dye]==false)
			{
				for (emn=0;emn<nemn;emn++) {
					f=0;
					for(exc=lowerexc[emn]; exc<upperexc[emn]; exc++){
						for (x=0;x<noPix;x++) {
							f+=Q[exc][dye]*Asub[x][dye]*Xsub[exc][x][emn]/QASsub[exc][x][emn];
						}}
					S[emn][dye]*=f/sumRowQA[emn][dye];	
					if (S[emn][dye]<spec_nothing) S[emn][dye]=spec_nothing;
				}
			}
		}
	}

	/**
	 * Do the update for the excitation spectra
	 */
	private void updateQ_NTF()
	{
		int x; int emn; int dye, exc;
		float[][]sumRowAS = new float[nexc][ndyes];
		float f;
		// calc estimated signal Q*A*S for the subspectra
		calc_QASsub();

		//calculate normalizer for the update rule
		//again, the normalizer is degenerate since it doesn't depend on exc
		//the implicit dependence on exc is due to the fact that some channels might not be 
		//illuminated in every excitation
		for (dye=0;dye<ndyes;dye++){
			for (exc=0;exc<nexc;exc++) {
				sumRowAS[exc][dye]=0;
				for (emn=0; emn<nchannels[exc]; emn++){
					for (x=0;x<noPix;x++){
						sumRowAS[exc][dye]+=S[channeltoemission[exc][emn]][dye]*Asub[x][dye];
					}}}}
		//do update
		for (dye=0;dye<ndyes;dye++) {
				for (exc=0;exc<nexc;exc++) {
					f=0;
					for(emn=0; emn<nchannels[exc]; emn++){
						for (x=0;x<noPix;x++) {
							f+=S[channeltoemission[exc][emn]][dye]*Asub[x][dye]*Xsub[exc][x][channeltoemission[exc][emn]]/QASsub[exc][x][channeltoemission[exc][emn]];
						}}
					Q[exc][dye]*=f/sumRowAS[exc][dye];				
			}
		}
	}

	
	/**
	 * Do the update for the concentrations
	 */
	private void updateA_NTF()
	{
		int x; int emn; int dye, exc;
		double f;
		float[] sumColQS = new float[ndyes];
		double onenorm, twonorm;

		// calc estimated signal A*S
		calc_QASsub();
		
		if (segbias>0)		//calculate the segregation bias term only if needed
		{
			for (x=0;x<noPix;x++)
			{
				for (dye=0, onenorm=0, twonorm=0;dye<ndyes;dye++)
				{
					f=Asub[x][dye];
					onenorm+=f;				// sum x_i
					twonorm+=f*f;			// sum x_i^2
				}
				onenorm/=twonorm;			
				twonorm=Math.sqrt(twonorm);	//twonorm is now the proper two norm sqrt(sum(x_i^2))
				onenorm/=twonorm;			//onenorm is now sum(x_i)/(sum x_i^2)^(3/2)
				for (dye=0;dye<ndyes;dye++)
				{							//calculate the gradient of |x|/||x|| and multiply by segbias
					Segbiasterm[x][dye]=segbias*(onenorm*Asub[x][dye]-1.0/twonorm);
				}		
			}
		}

		//calculate normalizer for update rule
		for (dye=0;dye<ndyes;dye++){
			sumColQS[dye]=0;
			for(exc=0; exc<nexc; exc++){
				for (emn=0;emn<nchannels[exc];emn++){
					sumColQS[dye]+=Q[exc][dye]*S[channeltoemission[exc][emn]][dye];
				}}}
		//do update
		for (dye=0;dye<ndyes;dye++) {
			for (x=0;x<noPix;x++) {
				f=0;
				for(exc=0; exc<nexc; exc++){
				for (emn=0;emn<nchannels[exc];emn++) {
					f+=Xsub[exc][x][channeltoemission[exc][emn]]/QASsub[exc][x][channeltoemission[exc][emn]]*S[channeltoemission[exc][emn]][dye]*Q[exc][dye];
				}}
				Asub[x][dye] *= (f+Segbiasterm[x][dye])/sumColQS[dye];
				if (Asub[x][dye]<conc_nothing) Asub[x][dye]=conc_nothing;	//make sure nothing becomes negative
			}
		}
	}


	/**
	 * Normalize the concentration to unity to pick up the differences in excitation efficiencies among spectra
	 * @param exc
	 */
	private void normalizeNMF(int exc){
		double mean;
		int emn, pix, dye;
		for(dye=0; dye<ndyes; dye++){
			mean=0;
			for(pix=0; pix<noPix; pix++){
				mean+=Asub[pix][dye];
			}
			mean*=0.01/noPix;
			for(pix=0; pix<noPix; pix++){
				Asub[pix][dye]/=mean;
			}
			for(emn=0; emn<nchannels[exc]; emn++){
				Sexc[exc][emn][dye]*=mean;
			}
		}
	}

	/**
	 * Normalize the concentration to unity to pick up the differences in excitation efficiencies among spectra
	 * @param exc
	 */
	private void normalizeNTF(){
		int dye, pix;
		int emn, exc;
		double temp;
		for (dye=0;dye<ndyes;dye++) {				//for each dye, calculate the 'E'-norm of the spectrum
			for (emn=0,power[dye]=0;emn<nemn;emn++)
				power[dye]+=Math.pow(S[emn][dye],E);
			power[dye]=(float) Math.pow(power[dye], 1f/E);
			for (emn=0;emn<nemn;emn++)
				S[emn][dye]/=power[dye];
		}
		if (E==1)  emnNorm="Normalized to unit area.";
		else emnNorm="Normalized to "+(new Float(E).toString())+" norm =1";
		
		for (dye=0;dye<ndyes;dye++) {				//for each dye, calculate the 'E'-norm of the spectrum
			temp=0;
			if (excNormType=="Max"){
				for (exc=0;exc<nexc;exc++)
					if (temp<Q[exc][dye]) temp=Q[exc][dye];
			}else{
				for (exc=0;exc<nexc;exc++)
					temp+=Math.pow(Q[exc][dye],E);
			}
			power[dye]*=(float) Math.pow(temp, 1f/E);
			for (exc=0;exc<nexc;exc++)
				Q[exc][dye]/= Math.pow(temp, 1f/E);			
		}
		if (excNormType=="Max")  excNorm="Normalized to maximum=1.";
		else excNorm="Normalized to "+(new Float(E).toString())+" norm=1";

		
		for (dye=0;dye<ndyes;dye++) {			//rescale S and A accordingly
			for (pix=0;pix<noPix;pix++)
				Asub[pix][dye]*=power[dye];
		}
		concNorm="Unnormalized.";
	}

	/**
	 * Normalize the concentration to unity to pick up the differences in excitation efficiencies among spectra
	 * @param exc
	 */
	private void normalizeNTF_maxConc(){
		int dye, pix;
		int emn, exc;
		double temp;
		for (dye=0;dye<ndyes;dye++) {				//for each dye, calculate the 'E'-norm of the spectrum
			for (emn=0,power[dye]=0;emn<nemn;emn++)
				power[dye]+=Math.pow(S[emn][dye],E);
			power[dye]=(float) Math.pow(power[dye], 1f/E);
			for (emn=0;emn<nemn;emn++)
				S[emn][dye]/=power[dye];
		}
		if (E==1)  emnNorm="Normalized to unit area.";
		else emnNorm="Normalized to "+(new Float(E).toString())+" norm =1";

		for (dye=0;dye<ndyes;dye++) {			//rescale S and A accordingly
			temp=0;
			for (pix=0;pix<noPix;pix++)
				if (Asub[pix][dye]>temp) temp=Asub[pix][dye];
			power[dye]*=temp;
			for (pix=0;pix<noPix;pix++)
				Asub[pix][dye]/=temp;
			
		}
		concNorm="Normalized to maximum approx 1";

		for (dye=0;dye<ndyes;dye++) {				//for each dye, calculate the 'E'-norm of the spectrum
			for (exc=0;exc<nexc;exc++)
				Q[exc][dye]*=power[dye];
		}
		excNorm="Unnormalized.";
	}

	/**
	 * init Asub with least square solution using current spectra
	 */
	private void solveforA() {
		int x;int emn, exc, dye, channel;
		calc_pinv();
		for (x=0;x<noPix;x++)		//loop over pixels and assign concentrations with random numbers
		{
			for (dye=0; dye<ndyes; dye++)
			{
				Asub[x][dye]=0;
				channel=0;
				for(exc=0; exc<nexc; exc++){
					for (emn=0; emn<nchannels[exc]; emn++){
						Asub[x][dye]+=pinvS[dye][channel]*Xsub[exc][x][channeltoemission[exc][emn]];
						channel++;
					}
				}
			}
		}
	}

	
	/**
	 * Prompt user for input, get background and start spectra
	 */
	private int parameterDialog()
	{
		GenericDialog nsources= new GenericDialog("Poisson NTF");
		nsources.addNumericField("Number of Sources", ndyes, 0);
		nsources.showDialog();
		if (nsources.wasCanceled())
		{	
			return CANCELLED;
		}
		else
		{
			ndyes=(int)nsources.getNextNumber();
			if (ndyes<1 || ndyes>10) 
			{
				IJ.error("Bad number of Sources!");
				return -1;
			}
			//
			allocate_ndyes_dependent();
			//get preferences from IJ preferences file
			getPrefs();

			GenericDialog pntf_dialog=new GenericDialog("Poisson NTF");
			pntf_dialog.addNumericField("_Number of Iterations", maxit, 0);
			pntf_dialog.addNumericField("Subsamples", subsamples, 0);
			pntf_dialog.addNumericField("Segregation Bias", segbias, 0);
			pntf_dialog.addNumericField("Saturation Threshold", saturation_threshold, 0);
			pntf_dialog.addNumericField("_Background Threshold", bg_threshold, 0);
			add_spectra_choice(pntf_dialog);
			pntf_dialog.addMessage("\n");
			pntf_dialog.addCheckbox("Specify Spectral Channels?", false);
			pntf_dialog.setOKLabel("Run!");

			pntf_dialog.showDialog();
			if (pntf_dialog.wasCanceled())
			{	
				return CANCELLED;
			}
			else
			{
				//RUNTIME
				//number of iterations
				maxit=(int) pntf_dialog.getNextNumber();
				if (maxit<1 || maxit>10000) IJ.error("Bad number of Iterations!");
				//Subsampling
				subsamples=(int) pntf_dialog.getNextNumber();
				if (subsamples<1) subsamples=1;
				//Segregation bias
				segbias=(float) pntf_dialog.getNextNumber();

				//INITIAL CONDITIONS
				//Saturation Threshold
				saturation_threshold=(float) pntf_dialog.getNextNumber();
				if (saturation_threshold<0 ) IJ.error("Threshold has to be positive!");
				//Background Threshold
				bg_threshold=(float) pntf_dialog.getNextNumber();
				if (bg_threshold<0 || bg_threshold>saturation_threshold) IJ.error("Lower threshold has to be positive\n" +
				"and below saturation!");
				determine_channel_order();
				gauss_spectra();
				read_spectra_choice(pntf_dialog);
				splice_spectra();
				//Spectral Windows
				if (pntf_dialog.getNextBoolean()) 
				{
					enter_emission_channels();
					determine_channel_order();
				}
			}
		}
		setPrefs();
		return 0;
	}

	/**
	 * init S with Gaussians
	 */
	private void gauss_spectra() {
		int emn,ii;
		int dye;
		//spectra - Gauss shape
		for (dye=0;dye<ndyes;dye++) {		//loop over dyes
			spectra_fixed[dye]=false;
			float sigma = (nemn-1f)/(ndyes+1f)+.1f;	//variance of gauss spectra
			for (emn=1;emn<=nemn;emn++) 				//loop over channels and assign spectra
			{
				ii=channel_order[emn-1];
				S[ii][dye]= (1f/(sigma*Math.sqrt(2*Math.PI))*Math.exp(-.5*Math.pow(emn-1-((dye+1)*nemn)/(ndyes+1f),2)/Math.pow(sigma,2))*(channel_width.get(ii)));
			}
		}
	}

	/**
	 * function that lets the user choose several stack, each of which is supposed to 
	 * correspond to one excitation
	 * @return
	 */
	
	private int openStacks()
	{
		//prompt user for the number excitation, i.e. the number of stacks to be opened
		GenericDialog nexcitations= new GenericDialog("Poisson NTF");
		nexcitations.addNumericField("Number of Excitations", ndyes, 0);
		nexcitations.showDialog();
		if (nexcitations.wasCanceled())
		{	
			return CANCELLED;
		}
		else
		{
			excitation_wavelength=new Vector<String>();
			nexc=(int) nexcitations.getNextNumber();
			allocate_nexc_dependent();
			if (nexc<1) {
				IJ.showMessage("Number of excitations has to be greater then 0!");
				return CANCELLED;
			}
			//loop over different excitations
			for (int exc=0; exc<nexc; exc++){
				IJ.run("Open...");
				imp[exc]=WindowManager.getCurrentImage();
				dim[exc]=imp[exc].getDimensions();
				//if stack is ok, convert to 32 gray scale and assign dimensions.
				if ((imp[exc].getBytesPerPixel()<32)) new StackConverter(imp[exc]).convertToGray32();
				excitation_wavelength.add(Integer.toString(exc));
				X[exc] = imp[exc].getStack();
				if (exc>0 && (w!=X[exc].getWidth() || h != X[exc].getHeight())) {
					IJ.showMessage("Stacks have to be of the same size!");
				}else{
					w = X[exc].getWidth();
					h = X[exc].getHeight();
					wh = w*h;	//number of pixels per frame
				}
			}
			//determine unique channels
			
			//make a list of the channel settings for each excitation
			for (int exc = 0; exc < nexc; exc++) {
				double []mid_wavelength=new double [X[exc].getSize()];
				if ((imp[exc].getOriginalFileInfo().fileName).endsWith(".lsm")){
					for (int ch = 0; ch < X[exc].getSize(); ch++) {
						mid_wavelength[ch]=Double.valueOf(X[exc].getSliceLabel(ch+1));
					}
				}else{
					for (int ch = 0; ch < X[exc].getSize(); ch++) {
						mid_wavelength[ch]=450+ch*50;
					}					
				}
				double channel_width;
				for (int ch = 0; ch < X[exc].getSize()-1; ch++) {
					channel_width=mid_wavelength[ch+1]-mid_wavelength[ch];
					ChannelsL.get(exc).add(mid_wavelength[ch]-0.5*channel_width);
					ChannelsR.get(exc).add(mid_wavelength[ch]+0.5*channel_width);
				}
				channel_width=mid_wavelength[X[exc].getSize()-1]-mid_wavelength[X[exc].getSize()-2];
				ChannelsL.get(exc).add(mid_wavelength[X[exc].getSize()-1]-0.5*channel_width);
				ChannelsR.get(exc).add(mid_wavelength[X[exc].getSize()-1]+0.5*channel_width);
			}

			//determine the number of unique emission channels.
			emission_wavelength=new Vector<Double>();
			channel_width=new Vector<Double>();
			double wl, width;
			for (int exc = 0; exc < ChannelsL.size(); exc++) {
				nchannels[exc]=ChannelsL.get(exc).size();
				for(int emn=0; emn < nchannels[exc]; emn++) {
					wl=0.5*(ChannelsL.get(exc).get(emn)+ChannelsR.get(exc).get(emn));
					width=(ChannelsR.get(exc).get(emn)-ChannelsL.get(exc).get(emn));
					if (!emission_wavelength.contains(wl)) {
						emission_wavelength.add(wl);
						channel_width.add(width);
					}
				}
			}
			nemn=emission_wavelength.size();
			//determined the number of emission channels, allocate memory
			allocate_nemn_dependent();
			for (int exc = 0; exc < ChannelsL.size(); exc++) {
				for(int emn=0; emn < nchannels[exc]; emn++) {
					wl=0.5*(ChannelsL.get(exc).get(emn)+ChannelsR.get(exc).get(emn));
					channeltoemission[exc][emn]=emission_wavelength.indexOf(wl);
				}
			}

			determine_exc_channel_map();
		}
		return 0;
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
		//get directory of experiment
		String LeicaDir = od.getDirectory();
		//get leica configuration file and read number of excitations and emissions from it
		File dir=new File(LeicaDir);
		String[] txt_list=dir.list(txt_file_filter);
		config_file=LeicaDir+txt_list[0];
		read_excitations_LeicaSP2(config_file);


		String series;
		try{
			series=LeicaFile.substring(LeicaFile.indexOf("Series"),LeicaFile.indexOf("Series")+9);
		}catch (Exception e)
		{
			IJ.showMessage("Unable to extract series from file name");
			return CANCELLED;
		}
		// read all images that of the same series and add them to a newly generated stack
		String[] im_list=dir.list(im_file_filter);
		int exc=0, channel=0, iminseries=0;
		for (int i=0; i<im_list.length; i++) {
			IJ.showProgress(i, im_list.length);
			if (im_list[i].contains(series)){
				if (imageselected[iminseries]==true) 
				{
					ImagePlus im=IJ.openImage(LeicaDir+im_list[i]);
					w=im.getWidth(); h=im.getHeight();
					if (X[exc]==null) X[exc]=new ImageStack(w,h);
					X[exc].addSlice("Channel "+(channel+1)+": "+ChannelsL.get(exc).get(channel)+'-'+ChannelsR.get(exc).get(channel)+"nm", im.getProcessor());
					channel++;
					if (channel==nchannels[exc]){exc++; channel=0;}
					if (exc==nexc) break;
				}
				iminseries++;
			}
		}
		//Display the stack and remember the data configuration file, assuming it is the only
		//text file in the directory
		for (exc=0; exc<nexc; exc++) {
			imp[exc] = new ImagePlus("Excitation "+(exc+1)+": "+excitation_wavelength.get(exc), X[exc]);
			imp[exc].show();
		}
		return 0;
	}

	/**
	 * calculate the center wavelength of all channels -- assumes the minimal and maximal
	 * wavelength entered by the user are center wavelength of first and last channel
	 * @param lambda_min
	 * @param lambda_max
	 * @return
	 */
	private int enter_emission_channels()
	{
		String msg="Enter emission channels.";
		GenericDialog enter_channels=new GenericDialog(msg);
		for (int emn=0; emn<nemn; emn++)
		{	//add numeric field with the current spectrum as default value
			enter_channels.addNumericField("Channel "+Integer.toString(emn+1)+" lower boundary: ", emission_wavelength.get(emn)-0.5*channel_width.get(emn), 6);
			enter_channels.addNumericField("Channel "+Integer.toString(emn+1)+" upper boundary: ", emission_wavelength.get(emn)+0.5*channel_width.get(emn), 6);
		}
		enter_channels.showDialog();
		if (enter_channels.wasCanceled())
		{	
			return -1;
		}
		else
		{
			double max, min;
			for (int emn=0; emn<nemn; emn++)
			{	//get user input
				min=enter_channels.getNextNumber();
				max=enter_channels.getNextNumber();
				emission_wavelength.set(emn, 0.5*(min+max));
				channel_width.set(emn, max-min);
			}
		}
		return 0;
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
		for (int dye = 0; dye < ndyes; dye++) {
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
		String [] labels=new String [ndyes];
		choose_spectra.addMessage("Keep spectra of dyes fixed?");
		for (int dye = 0; dye < ndyes; dye++) 
		{
			labels[dye]="Dye ";
			labels[dye]=labels[dye].concat(Integer.toString(dye+1));
		}
		choose_spectra.addCheckboxGroup(1, ndyes, labels, spectra_fixed);
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
		for (int dye = 0; dye < ndyes; dye++) {
			spec_choice[dye]=choose_spectra.getNextChoice();
			if (spec_choice[dye].equals(items[0]) || spec_choice[dye].equals(items[3])) {}	//Gaussian
			else if (spec_choice[dye].equals(items[1])) getSpectrum(dye);			//ROI
			else if (spec_choice[dye].equals(items[2])) enterSpectrum(dye);		//Manually
			else readSpectrum(dye, curDir, spec_choice[dye]);				//from file
		}
		//read whether spectra are to be kept fixed
		for (int dye = 0; dye < ndyes; dye++) {
			spectra_fixed[dye]=choose_spectra.getNextBoolean();
		}

		//Normalize the spectra
		for (int dye=0;dye<ndyes;dye++) {				//for each dye, calculate the 'E'-norm of the spectrum
			power[dye]=0.0;
			for (int emn=0;emn<nemn;emn++)
				power[dye]+=Math.pow(S[emn][dye],E);
			power[dye]=(float) Math.pow(power[dye], 1f/E);
		}
		for (int dye=0;dye<ndyes;dye++) {			//rescale S
			for (int emn=0;emn<nemn;emn++)
			{
				S[emn][dye]/=power[dye];
			}
		}
		//save initial spectra for nostalgia and normalize them for display
		int ii;
		for (int dye = 0; dye < ndyes; dye++) {
			for (int emn=0; emn<nemn; emn++)	
			{
				ii=inverse_channel_order[emn];
				initialS[emn][dye]=(float) S[ii][dye];
			}
		}
		float maxS=0;
		for (int dye=0;dye<ndyes;dye++) {
			for (int emn=0;emn<nemn;emn++) {
				ii=inverse_channel_order[emn];
				if (initialS[emn][dye]/channel_width.get(ii)>maxS) 
					maxS=(float) ((float) initialS[emn][dye]/channel_width.get(ii));
			}
		}
		for (int dye=0;dye<ndyes;dye++) {
			for (int emn=0;emn<nemn;emn++) {
				ii=inverse_channel_order[emn];
				initialS[emn][dye]=(float) (initialS[emn][dye]/channel_width.get(ii)/maxS);
			}
		}
		
		for(int dye=0; dye<ndyes; dye++)
		{
			for (int exc=0;exc<nexc;exc++){
				for (int emn=0;emn<nchannels[exc];emn++) 				//loop over channels and assign spectra
				{
					Sexc[exc][emn][dye]=S[channeltoemission[exc][emn]][dye];
				}
			}
		}
		return 0;
	}
	/**
	 * User interface to input background signal, either by ROI or by number 
	 */
	private int getBackground(String choice)
	{
		int emn, exc;
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

		int goodinput=0;
		if (choice.equals(bgitems[1]))		//open a dialog asking the user to select a ROI
		{
			new WaitForUserDialog("Background Selection", "Select Background ROI").show();
			Roi	bg_roi=null; 
			int currentslice_signal=0;
			for (exc=0; exc<nexc; exc++){
				if (imp[exc].getRoi()!=null) 
				{
					bg_roi=imp[exc].getRoi();
					//determine current slice of the image to extract the background
					currentslice_signal=imp[exc].getCurrentSlice();
				}}
			currentslice_signal--;	//correct for counting with base 1
			currentslice_signal-=currentslice_signal%nemn;	//determine first channel

			if (bg_roi==null)
				goodinput=NOROISELECTED;	//no ROI selected -- have the user try again
			else
				goodinput=GOODINPUT;

			if (goodinput==GOODINPUT)	//read ROI, calculate background intensity
			{
				for(exc=0; exc<nexc; exc++){
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
								for (emn=0; emn<nchannels[exc]; emn++)
								{
									v=((float [])X[exc].getPixels(currentslice_signal+emn+1))[y*w+x];
									bg[exc][emn]+=v;
									bg_sigma[exc][emn]+=v*v;	
								}
							}
						}
					}
					for (emn=0; emn<nchannels[exc]; emn++)
					{
						bg[exc][emn]/=(float)no_bgpix;
						bg_sigma[exc][emn]/=(float)no_bgpix;
						bg_sigma[exc][emn]=(float) Math.sqrt(bg_sigma[exc][emn]-bg[exc][emn]*bg[exc][emn]);
					}
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
				for (exc=0; exc<nexc; exc++)
					for (emn=0; emn<nemn; emn++)
					{
						bg[exc][emn]=bg_uniform; bg_sigma[exc][emn]=1;
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

	private void splice_spectra()
	{
		int emn1, emn2, exc1, exc2, dye;
		double [][][] excitation1=new double [nexc][nexc][ndyes];
		double [][][] excitation2=new double [nexc][nexc][ndyes];
		for (exc1 = 0; exc1 < nexc; exc1++) {
			for (exc2=0; exc2<exc1; exc2++){
				for (dye=0; dye<ndyes; dye++){
					excitation1[exc1][exc2][dye]=0;								
					excitation2[exc1][exc2][dye]=0;								
				}
			}
		}
		for (exc1 = 0; exc1 < nexc; exc1++) {
			for (exc2=0; exc2< nexc; exc2++){
				for (emn1=0; emn1<nchannels[exc1]; emn1++){
					for (emn2=0; emn2<nchannels[exc2]; emn2++){
						if (channeltoemission[exc1][emn1]==channeltoemission[exc2][emn2]){
							for (dye=0; dye<ndyes; dye++){
								excitation1[exc1][exc2][dye]+=Sexc[exc1][emn1][dye];								
								excitation2[exc1][exc2][dye]+=Sexc[exc2][emn2][dye];								
							}
						}
					}
				}
			}
		}
		for (exc1 = 0;  exc1<nexc; exc1++) {
			for(dye=0; dye<ndyes; dye++){
				if (exc1==0) Q[exc1][dye]=1;
				else{
					Q[exc1][dye]=Q[exc1-1][dye]*excitation2[exc1-1][exc1][dye]/excitation1[exc1-1][exc1][dye];
				}
			}			
		}
		double [][]norm= new double [nemn][ndyes];
		for (dye=0; dye<ndyes; dye++){
			for (emn1=0; emn1<nemn; emn1++){
				S[emn1][dye]=0;
				norm[emn1][dye]=0;
			}
		}
		for (exc1=0; exc1<nexc; exc1++){
			for (dye=0; dye<ndyes; dye++){
				for (emn1=0; emn1<nchannels[exc1]; emn1++){
					S[channeltoemission[exc1][emn1]][dye]+=Sexc[exc1][emn1][dye];
					norm[channeltoemission[exc1][emn1]][dye]+=Q[exc1][dye];
				}
			}
		}
		for (dye=0; dye<ndyes; dye++){
			for (emn1=0; emn1<nemn; emn1++){
				S[emn1][dye]/=norm[emn1][dye];
			}
		}
	}

	/**
	 * function that prompts the user for a ROI and calculates the mean spectrum of the 
	 * pixels in that ROI
	 * @param dye
	 */
	private int getSpectrum(int dye)
	{
		int exc, emn;
		boolean goodroi;
		String msg="Select a ROI for dye ";
		msg=msg.concat(Integer.toString(dye+1));
		new WaitForUserDialog("Start Spectra", msg).show();
		Roi spec_roi=null;
		int currentslice_signal=0;
		for(exc=0; exc<nexc; exc++) {
			if (imp[exc].getRoi()!=null)
			{
				spec_roi=imp[exc].getRoi();
				currentslice_signal=imp[exc].getCurrentSlice();

			}
		}
		currentslice_signal--;	//correct for base 1 counting
		currentslice_signal-=currentslice_signal%nemn;	//determine slice of first channel.
		if (spec_roi==null) goodroi=false;
		else goodroi=true;
		if (goodroi)
		{
			//determine time and depth slice currently on display
			int npix=0; float v;
			for (exc=0; exc<nexc; exc++){
				for (emn=0; emn<nemn; emn++)
				{
					Sexc[exc][emn][dye]=0;
				}
			}
			for (int x=0; x<w; x++) //loop over all pixels (terribly inefficient)
			{
				for (int y=0; y<h; y++)
				{
					if (spec_roi.contains(x,y))	//calc mean from pixels in ROI
					{	
						npix++;
						for (exc=0; exc<nexc; exc++){
							for (emn=0; emn<nchannels[exc]; emn++)
							{
								v=((float [])X[exc].getPixels(currentslice_signal+emn+1))[y*w+x];
								Sexc[exc][emn][dye]+=v;
							}
						}
					}
				}
			}

			for (exc=0; exc<nexc; exc++) {
				for (emn=0; emn<nemn; emn++) {
					Sexc[exc][emn][dye]/=(float)npix;	//normalize
				}
			}
			for (exc=0; exc<nexc; exc++) 
				for (emn=0; emn<nchannels[exc]; emn++)
					Sexc[exc][emn][dye]-=bg[exc][emn];		//subtract previously determined background
			return GOODINPUT;
		}
		else 
		{	//complain since user did not select a ROI
			IJ.showMessage("Select ROI!");
			return NOROISELECTED;
		}
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
			int emn=0;
			double step, result;	//integration step
			for (emn = 0; emn < nemn; emn++) 
			{
				step=channel_width.get(emn)/20;	//integration step
				result=0;
				for (double lambdax=emission_wavelength.get(emn)-0.5*channel_width.get(emn);lambdax<emission_wavelength.get(emn)+0.5*channel_width.get(emn); lambdax+=step)
				{	//use midpoint of interval as pivot for interpolated value of spectrum
					result+=interpolate(lambda, emission, lambdax+0.5*step)*step;
				}
				S[emn][dye]=result;
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
			//result=Math.min(F1, F2);
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
	 * loop over all slices and pixel and determine the minimal signal for each channel
	 */
	private void minimumBackground()
	{
		float min;
		int currentslice_signal;
		for (int exc=0; exc<nexc; exc++)
		{
			for (int emn=0; emn<nchannels[exc]; emn++)
			{
				for (int slice=0; slice<tz_dim; slice++)
				{
					currentslice_signal=slice*nemn+emn+1;
					min=(float) saturation_threshold;
					for (int pix=0; pix<wh; pix++)
					{	//determine minimum in each channel
						if (((float [])X[exc].getPixels(currentslice_signal))[pix]<min) min=((float [])X[exc].getPixels(currentslice_signal))[pix];
					}
					bg[exc][emn]=min;
					bg_sigma[exc][emn]=0;
				}
			}
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
		msg="Enter spectrum of dye ";
		msg=msg.concat(Integer.toString(dye+1));
		GenericDialog enter_spectra=new GenericDialog(msg);
		for (int emn=0; emn<nemn; emn++)
		{	//add numeric field with the current spectrum as default value
			enter_spectra.addNumericField(Integer.toString(emn+1), S[emn][dye], 4);
		}
		enter_spectra.showDialog();
		if (enter_spectra.wasCanceled())
		{	
			return CANCELLED;
		}
		else
		{
			for (int emn=0; emn<nemn; emn++)
			{	//get user input
				S[emn][dye]=enter_spectra.getNextNumber();
				if (S[emn][dye]<spec_nothing) S[emn][dye]=spec_nothing;
			}
		}
		for (int exc=0;exc<nexc;exc++){
			for (int emn=0;emn<nchannels[exc];emn++) 				//loop over channels and assign spectra
			{
				Sexc[exc][emn][dye]=S[channeltoemission[exc][emn]][dye];
			}
		}
		return GOODINPUT;
	}

	/**
	 * decide which pixels are saturated or background
	 */
	private void calcBgMask()
	{
		int currentpix, currentslice_signal;
		int emn,exc;
		for (int slice=0; slice<tz_dim; slice++)
		{
			currentpix=slice*wh;
			//initialize the background map with n
			for (int pix=0; pix<wh; pix++)
				bgMask[pix+currentpix]=totalchannels;

			//background code: each saturated channel subtracts 2n from bgMask -> a single channel makes it negative
			//each channel below background[lambda]+bg_threshold subtracts 1 from bgMask -> it remains positive unless
			//all channels are below bg
			for (exc=0; exc<nexc; exc++)
			{
				for (emn=0; emn<nchannels[exc]; emn++)
				{
					currentslice_signal=slice*nemn+emn+1;
					for (int pix=0; pix<wh; pix++)
					{	//for each pixel, decide whether it is below background or above saturation
						if (((float [])X[exc].getPixels(currentslice_signal))[pix]<bg[exc][emn]+no_std*bg_sigma[exc][emn]+bg_threshold) bgMask[currentpix+pix]--;
						if (((float [])X[exc].getPixels(currentslice_signal))[pix]>saturation_threshold) bgMask[currentpix+pix]-=2*totalchannels;
					}
				}
			}
		}
	}


	private void determine_exc_channel_map(){
		int exc, emn, emn1;
		int min, max;
		for (emn=0; emn<nemn; emn++){
			min=nexc; max=0;
			for (exc=0; exc<nexc; exc++){
				for (emn1=0; emn1<nchannels[exc]; emn1++){
					if (channeltoemission[exc][emn1]==emn){
						if (exc<min) min=exc;
						if (exc>max) max=exc;
					}
				}
			}
			lowerexc[emn]=min;
			upperexc[emn]=max+1;
		}
	}

	/**
	 * Function that parses the Leica SP2 configuration file
	 * @param conf_file
	 */
	private int read_excitations_LeicaSP2(String conf_file)
	{
		int NPMTS=4;
		int exc, emn, i;
		excitation_wavelength=new Vector<String>();
		try {
			BufferedReader LeicaInfoFile=new BufferedReader(new FileReader(conf_file));
			String line;
			StringTokenizer entries;
			Vector <String> temp_excitation_wavelength=new Vector<String>();
			Vector <Integer> aquisition=new Vector <Integer>();			//vector to hold the excitation to which the acquisition sequence corresponds to
			Vector<Vector<Double>> SPL=new Vector<Vector<Double>>();	//save the settings for each mirror and image sequence
			Vector<Vector<Double>> SPR=new Vector<Vector<Double>>();
			for (i = 0; i < NPMTS; i++) {
				SPL.add(null); SPL.set(SPL.size()-1, new Vector<Double>());
				SPR.add(null); SPR.set(SPR.size()-1, new Vector<Double>());
			}
			//read in all line from the file and assume the first column is wavelength
			//the second column is spectrum
			int sequence=1;
			String seq_base="HARDWARE PARAMETER #0 SEQUENCE ";
			//read everything in the file until the parameters for the first data sequence are found 
			while((line=LeicaInfoFile.readLine())!=null &&!(line).contains(seq_base+sequence))
			{}
			temp_excitation_wavelength.add("");
			int channel=0;
			while((line=LeicaInfoFile.readLine())!=null)
			{
				if (line.contains("AOTF"))
				{
					entries=new StringTokenizer(line);
					entries.nextToken();
					String wavelength=entries.nextToken();
					String value=entries.nextToken();
					double aotfval=0; 
					//read entry where a number is expected
					try{aotfval=Double.valueOf(value);}
					catch (Exception e) {;}
					if (aotfval>0.01) 
					{
						temp_excitation_wavelength.set(sequence-1, temp_excitation_wavelength.get(sequence-1)+wavelength);
					}
				}
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
					if (lr.equals("(left)")) SPL.get(channel).add(wl);
					else if (lr.equals("(right)")) {SPR.get(channel).add(wl); channel++;}
				}
				if (line.contains(seq_base)) 
				{sequence++; temp_excitation_wavelength.add(""); channel=0;	}
			}
			
			//TODO:
			int detectedexc=temp_excitation_wavelength.size();
			GenericDialog selectimages= new GenericDialog("Poisson NTF");
			imageselected=new boolean [NPMTS*SPL.get(0).size()+1];
			String [] labels=new String [detectedexc];
			for (int npmt=0; npmt<NPMTS; npmt++)
			{
				for (exc=0; exc<SPL.get(npmt).size(); exc++)
				{
					imageselected[NPMTS*exc+npmt]=true;
				}
			}
			for (exc=0; exc<detectedexc; exc++)
			{
				labels[exc]=temp_excitation_wavelength.get(exc);
			}
			selectimages.addCheckboxGroup(detectedexc, 1, labels, imageselected);
			selectimages.showDialog();
			if (selectimages.wasCanceled())
			{	
				return CANCELLED;
			}else{
				int alreadyremoved=0;
				for (exc=0; exc<detectedexc; exc++)
				{	
					if (selectimages.getNextBoolean()==false)
					{
						temp_excitation_wavelength.remove(exc-alreadyremoved);
						sequence--;
						for (int npmt=0; npmt<NPMTS; npmt++)
						{
							imageselected[NPMTS*exc+npmt]=false;
							SPL.get(npmt).remove(exc-alreadyremoved);
							SPR.get(npmt).remove(exc-alreadyremoved);
						}
						alreadyremoved++;
					}
				}
			}
			
			//determine the number of unique excitations. 
			nexc=0;
			for (exc=0; exc<temp_excitation_wavelength.size(); exc++)
			{
				int unique=0, first=-1;
				for (int exc1=0; exc1<exc; exc1++)
				{
					if (temp_excitation_wavelength.get(exc).equals(temp_excitation_wavelength.get(exc1))) 
					{
						unique++;
						if (first<0) first=exc1;
					}
				}
				//if (unique==0) {excitation_wavelength.add(temp_excitation_wavelength.get(exc)); aquisition.add(exc); nexc++;}
				//else aquisition.add(first);
				excitation_wavelength.add(temp_excitation_wavelength.get(exc)); aquisition.add(exc); nexc++;
			}
			//determined number of excitations, allocate stuff to do with excitations.
			allocate_nexc_dependent();

			//make a list of the channel settings for each excitation
			for (int acq = 0; acq < sequence; acq++) {
				for (int pmt = 0; pmt < NPMTS; pmt++) {
					ChannelsL.get(aquisition.get(acq)).add(SPL.get(pmt).get(acq));
					ChannelsR.get(aquisition.get(acq)).add(SPR.get(pmt).get(acq));
				}
			}
			//determine the number of unique emission channels.
			emission_wavelength=new Vector<Double>();
			channel_width=new Vector<Double>();
			double wl, width;
			for (exc = 0; exc < ChannelsL.size(); exc++) {
				nchannels[exc]=ChannelsL.get(exc).size();
				for(emn=0; emn < nchannels[exc]; emn++) {
					wl=0.5*(ChannelsL.get(exc).get(emn)+ChannelsR.get(exc).get(emn));
					width=(ChannelsR.get(exc).get(emn)-ChannelsL.get(exc).get(emn));
					if (!emission_wavelength.contains(wl)) {
						emission_wavelength.add(wl);
						channel_width.add(width);
					}
				}
			}
			nemn=emission_wavelength.size();
			//determined the number of emission channels, allocate memory
			allocate_nemn_dependent();
			for (exc = 0; exc < ChannelsL.size(); exc++) {
				for(emn=0; emn < nchannels[exc]; emn++) {
					wl=0.5*(ChannelsL.get(exc).get(emn)+ChannelsR.get(exc).get(emn));
					channeltoemission[exc][emn]=emission_wavelength.indexOf(wl);
				}
			}
			determine_exc_channel_map();

		} catch (FileNotFoundException e) {
			return 0;
		} catch (IOException e) {
			return 0;
		}
		return 0;
	}
	/**
	 * function that determines 
	 * the order of the channels with increasing wavelength
	 * this is important to produce reasonable spectra plots if the images in the source stack are
	 * not ordered, as it can happen with Leica SP2 data.
	 */
	private void determine_channel_order()
	{
		double prev_min=-1, min;
		int c=0;
		for (int emn=0; emn<nemn; emn++) 
		{
			min=1000;
			for(int emn1=0; emn1<nemn; emn1++)
			{
				if (emission_wavelength.get(emn1)<min && emission_wavelength.get(emn1)>prev_min)
				{min=emission_wavelength.get(emn1); c=emn1;}
			}
			channel_order[emn]=c;
			prev_min=min;
		}
		for (int emn=0; emn<nemn; emn++) 
		{
			inverse_channel_order[channel_order[emn]]=emn;
		}
	}

	/**
	 * generate new plot for the spectra
	 * @return 
	 * 
	 */
	public void plotEmissionSpectra(boolean newWindow) {
		int dye, emn,ii;
		double [][] spec=new double [nemn][ndyes];
		for (dye=0;dye<ndyes;dye++) {
			for (emn=0;emn<nemn;emn++) {
				ii=channel_order[emn];
				spec[emn][dye]=S[ii][dye];
			}
		}
		//colors for the different spectra
		Color[] colors = {Color.blue,Color.green, Color.red, Color.cyan, Color.gray, Color.darkGray};
		Plot plot = new Plot("PoissonNTF emission spectra. "+emnNorm,"wave length [nm]","intensity",(float[])null,(float[])null, PlotWindow.LINE);
		if (plotEmission!=null && newWindow){
			plotEmission.close();
		}

		if (plotEmission==null || newWindow){
			plotEmission=new PlotWindow("PoissonNTF emission spectra. "+emnNorm,"wave length [nm]","intensity",(float[])null,(float[])null);
		}
		float[] rowS = new float[nemn];
		float[] lambdas=new float [nemn];
		for (emn = 0; emn < lambdas.length; emn++) {
			ii=channel_order[emn];
			lambdas[emn]=(float) (double) (emission_wavelength.get(ii)); //channel wavelength as x-axis
		}
		plot.setLimits(emission_wavelength.get(channel_order[0]), emission_wavelength.get(channel_order[nemn-1]), 0, 1);
		int c=0; float maxS=0;
		//determine maximum of all spectra
		for (dye=0;dye<ndyes;dye++) {
			for (emn=0;emn<nemn;emn++) {
				ii=inverse_channel_order[emn];
				if (spec[emn][dye]/channel_width.get(ii)>maxS) 
					maxS=(float) ((float) spec[emn][dye]/channel_width.get(ii));
			}
		}
		//loop over spectra and plot each
		for (dye=0;dye<ndyes;dye++) {
			for (emn=0;emn<nemn;emn++) {
				ii=inverse_channel_order[emn];
				rowS[emn]=(float) (spec[emn][dye]/channel_width.get(ii)/maxS);
			}
			plot.setLineWidth(2);            
			plot.setColor(colors[c]);
			plot.addPoints(lambdas, rowS, PlotWindow.LINE);
			plot.addLabel(0.05, 0.1*(dye+1), "Dye : "+(new Float(dye+1)).toString());
			//plot initial spectra only if not kept fixed
			if (true) //(spectra_fixed[dye]==false)
			{
				for (emn=0;emn<nemn;emn++) {
					ii=inverse_channel_order[emn];
					rowS[emn]=(float) initialS[emn][dye];
				}
				plot.setLineWidth(1);            
				plot.addPoints(lambdas, rowS, PlotWindow.LINE);
			}
			if (c++>colors.length) c=0;
		}
		plot.draw();
		plotEmission.drawPlot(plot);
	}

	/**
	 * generate new plot for the spectra
	 * @return 
	 * 
	 */
	public void plotExcitationSpectra(boolean newWindow) {
		int dye, exc;
		double [][] spec=new double [nexc][ndyes];
		for (dye=0;dye<ndyes;dye++) {
			for (exc=0;exc<nexc;exc++) {
				spec[exc][dye]=Q[exc][dye];
			}
		}
		//colors for the different spectra
		Color[] colors = {Color.blue,Color.green, Color.red, Color.cyan, Color.gray, Color.darkGray};
		Plot plot = new Plot("PoissonNMF excitation spectra. "+excNorm,"Excitation","intensity",(float[])null,(float[])null, PlotWindow.LINE);
		if (plotExcitation!=null && newWindow){
			plotExcitation.close();
		}
		
		if (plotExcitation==null || newWindow){
			plotExcitation=new PlotWindow("PoissonNTF excitation spectra. "+excNorm,"Excitation","intensity",(float[])null,(float[])null);
		}
		float[] rowS = new float[nexc];
		float[] lambdas=new float [nexc];
		for (exc = 0; exc < lambdas.length; exc++) {
			lambdas[exc]=(float) exc+1;//channel wavelength as x-axis
		}
		int c=0; float maxS=0;
		//determine maximum of all spectra
		for (dye=0;dye<ndyes;dye++) {
			for (exc=0;exc<nexc;exc++) {
				if (spec[exc][dye]>maxS) 
					maxS=(float) spec[exc][dye];
			}
		}
		plot.setLimits(lambdas[0]-0.3, lambdas[exc-1]+0.3, 0, maxS);
		//loop over spectra and plot each
		for (dye=0;dye<ndyes;dye++) {
			for (exc=0;exc<nexc;exc++) {
				rowS[exc]=(float) spec[exc][dye];
			}
			plot.setLineWidth(2);            
			plot.setColor(colors[c]);
			plot.addLabel(0.05, 0.1*(dye+1), "Dye : "+(new Float(dye+1)).toString());
			plot.addPoints(lambdas, rowS, PlotWindow.LINE);
			if (c++>colors.length) c=0;
		}
		plot.draw();
		plotExcitation.drawPlot(plot);
	}

	/**
	 * Plot background spectrum
	 * @return
	 */
	public PlotWindow plotbg() {
		int emn, exc;
		//color fot the spectrum
		Color[] colors = {Color.blue, Color.green, Color.red};
		Plot plot=new Plot("Poisson NMF background spectra","wave length [nm]","intensity",(float[])null,(float[])null,PlotWindow.LINE);
		double[] lambdas=new double [nemn];
		for (emn = 0; emn < lambdas.length; emn++) {
			lambdas[emn]=emission_wavelength.get(emn); //channel wavelength as x-axis
		}
		double minbg=saturation_threshold;
		double maxbg=0;
		for (exc=0; exc<nexc; exc++)
		{
			for (emn = 0; emn < nchannels[exc]; emn++) {
				if (bg[exc][emn]>maxbg) maxbg=bg[exc][emn]; 
				if (bg[exc][emn]<minbg) minbg=bg[exc][emn];
			}
		}
		if (minbg==maxbg) {minbg-=5; maxbg+=5;}
		plot.setLimits(emission_wavelength.get(0), emission_wavelength.get(nemn-1), minbg, maxbg);

		//loop over spectra and plot each

		for (exc=0; exc<nexc; exc++)
		{
			float [] lambdach=new float [nchannels[exc]];
			float [] bgplot=new float [nchannels[exc]];
			for (emn = 0; emn < nchannels[exc]; emn++) {
				lambdach[emn]=(float) ((double) emission_wavelength.get(channeltoemission[exc][emn])); //channel wavelength as x-axis				
				bgplot[emn]=(float) bg[exc][emn];
			}
			plot.setLineWidth(2);            
			plot.setColor(colors[exc%3]);
			plot.addPoints(lambdach, bgplot, PlotWindow.LINE);
		}

		return plot.show();
	}

	
	/**
	 * show the RGB overlay of specified sources
	 */
	public void RGB_overlay()
	{
		if (currentnormalization==1){
			togglenormalization();
		}
		
		int [] sources=new int [3];
		int currentslice_conc;
		//prompt the user for an assignment of channels to Red, Green, Blue
		GenericDialog sourcesDialog=new GenericDialog("Select Sources");
		if (ndyes>2)
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
				currentslice_conc-=currentslice_conc%ndyes;	//determine the first channel
			}
			else currentslice_conc=0;
			sources[0]=(int) sourcesDialog.getNextNumber();
			sources[1]=(int) sourcesDialog.getNextNumber();
			if (ndyes>2) sources[2]=(int) sourcesDialog.getNextNumber();
			//check that entered source values are positive and not larger than the number of sources
			for (int i=0; i<Math.min(ndyes,3);i++)
				if (sources[i]>ndyes  || sources[i]<1) 
				{
					IJ.showMessage("Invalid source number!");
					return;
				}
			ImagePlus overlay=null;
			overlay=NewImage.createRGBImage("Overlay of concentrations", w,h,1, NewImage.FILL_BLACK);
			int source;
			//normalize images to max==1
			//produce the RGB overlay if the number of sources == 3
			if (ndyes>2)
			{
				for (int pix=0; pix<wh; pix++)
				{
					for (source=0; source<3; source++)
						((int[])overlay.getProcessor().getPixels())[pix]+=
							((int)(((1<<8)-1.0)*Math.min(Math.abs(((float[])A.getPixels(currentslice_conc+sources[source]))[pix]), 1.0)))<<(8*source);
				}
			}
			else
			{
				for (int pix=0; pix<wh; pix++)
				{
					for (source=0; source<2; source++)
						((int[])overlay.getProcessor().getPixels())[pix]+=
							((int)(((1<<8)-1.0)*Math.min(Math.abs(((float[])A.getPixels(currentslice_conc+sources[source]))[pix]), 1.0)))<<(8*(source+1));
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
			currentpix-=currentpix%ndyes;	//determine first color channel slice
			currentpix/=ndyes;				//number of current tz-slice
			currentpix*=wh;				//pixel offset = #tz-slices * pixel per frame
		}
		else currentpix=0;
		//construct the map of background, signal and saturated regions
		for (int pix=0; pix<wh; pix++)
		{
			if (bgMask[currentpix+pix]==2*totalchannels) 	//forground pixels
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
	 * plot/show image stack of concentrations (A)
	 * 
	 */
	private void showConcentrations(boolean doPostIter) {
		int dye, emn , exc, channel, count=0;

		segbias=0;
		if (doPostIter) for(int i=0; i<no_postiterations; i++) updateA_NTF();
		if (img!=null){
			img.changes=false;
			img.close();
		}
		//allocate a new image stack, one slice for each dye at each t and z slice
		A=new ImageStack(w,h);
		for (int j=0; j<tz_dim; j++)	//time and depth slices
			for (dye=0; dye<ndyes; dye++)		//add slice for each dye
				A.addSlice("source "+(dye+1), new FloatProcessor(w,h) );

		if (currentnormalization==1){
			normalizeNTF();
		}else{
			normalizeNTF_maxConc();
		}
		calc_pinv();
		double [] c=new double [ndyes];
		int currentslice_conc,currentslice_signal,  currentpix;
		//loop over all slices and pixels and put together the image stack of the NMF sources
		for (int slice=0; slice<tz_dim; slice++)
		{
			currentpix=slice*wh;
			currentslice_conc=slice*ndyes;

			for (int pix=0; pix<wh; pix++)
			{
				if (bgMask[currentpix+pix]==2*totalchannels) 	//fill in NMF pixels
				{
					for (dye=0; dye<ndyes; dye++)
						((float[])A.getPixels(currentslice_conc+dye+1))[pix]=(float)Asub[count][dye];
					count++;
				}
				else	//solve least squares for all other pixels
				{
					//construct vector, set negative values to signal_nothing
					//solve for concentrations
					for (dye=0; dye<ndyes; dye++)
					{
						c[dye]=0;
						channel=0;
						for(exc=0; exc<nexc; exc++){
							currentslice_signal=slice*nchannels[exc];
							for (emn=0; emn<nchannels[exc]; emn++){	
								c[dye]+=pinvS[dye][channel]*(((float[])X[exc].getPixels(currentslice_signal+emn+1))[pix]-bg[exc][emn]);
								channel++;
							}}
					}
					for (dye=0; dye<ndyes; dye++)
						((float[])A.getPixels(currentslice_conc+dye+1))[pix]=(float) c[dye];
				}
			}
		}
		//Display
		img = new ImagePlus("NTF sources. "+concNorm, A);
		if (isHyperstack) //for hyperstacks, set the time and depth dimensions
		{
			img.setDimensions(ndyes,dim[0][3], dim[0][4]);
			img.setOpenAsHyperStack(true);
		}
		img.show();
		img.updateAndDraw();
	}
	
	public void togglenormalization(){
		if (currentnormalization==0){
			currentnormalization=1;
		}else{
			currentnormalization=0;			
		}
		showConcentrations(false);
		plotEmissionSpectra(true);
		plotExcitationSpectra(true);
	}

	private void allocate_nexc_dependent()
	{
		//Vectors that save the emission channel boundaries for each excitation.
		ChannelsL=new Vector<Vector<Double>>();
		ChannelsR=new Vector<Vector<Double>>();
		for (int exc=0; exc<nexc; exc++)
		{
			ChannelsL.add(null); ChannelsL.set(exc, new Vector<Double>());
			ChannelsR.add(null); ChannelsR.set(exc, new Vector<Double>());
		}
		//dimension for each excitation
		dim =new int [nexc][5];
		nchannels=new int [nexc];
		//allocate image stacks for different excitations
		imp=new ImagePlus [nexc];
		X=new ImageStack [nexc];
	}

	private void allocate_ndyes_dependent() {
		S = new double [nemn][ndyes];
		Sexc = new double [nexc][nemn][ndyes];
		Q = new double [nexc][ndyes];
		initialS = new float [nemn][ndyes];
		spectra_fixed=new boolean [ndyes];		
		power=new double [ndyes];
		spec_choice=new String [ndyes];
		for(int exc=0; exc<nexc; exc++) totalchannels+=nchannels[exc];
		pinvS=new double [ndyes][totalchannels];
	}


	private void allocate_nemn_dependent()
	{
		channeltoemission= new int [nexc][nemn];
		bg=new double [nexc][nemn];
		bg_sigma=new double [nexc][nemn];
		channel_order=new int [nemn];
		inverse_channel_order=new int [nemn];
		lowerexc=new int [nemn];
		upperexc=new int [nemn];		
	}

	/**
	 * Function that reads the previous settings from the imageJ preferences file
	 **/
	private void getPrefs() {
		String label;
		maxit=(int)Prefs.get("PoissonNTF.maxit", maxit);
		segbias=Prefs.get("PoissonNTF.segbias", segbias);
		saturation_threshold=Prefs.get("PoissonNTF.saturation_threshold", saturation_threshold);
		bg_threshold=Prefs.get("PoissonNTF.bg_threshold", bg_threshold);
		bg_choice=Prefs.get("PoissonNTF.bg_choice", "none");
		subsamples=(int)Prefs.get("PoissonNTF.subsamples", subsamples);

		//choice of initial spectra and decision to keep some spectra fixed
		if ((int)Prefs.get("PoissonNTF.ndyes",0)==ndyes)
		{
			for(int dye=0; dye<ndyes; dye++)
			{			
				label="PoissonNTF.Dye_";
				label=label.concat(Integer.toString(dye+1));
				spec_choice[dye]=Prefs.get(label, "none");
				label="PoissonNTF.DyeFixed_";
				label=label.concat(Integer.toString(dye+1));
				spectra_fixed[dye]=Prefs.get(label, false);
			}
		}
	}

	/**
	 * Save the parameters to the imageJ preferences file: one-to-one correspondence
	 * to getPrefs() above
	 */
	private void setPrefs() {
		String label;
		Prefs.set("PoissonNTF.maxit", maxit);
		Prefs.set("PoissonNTF.segbias", segbias);
		Prefs.set("PoissonNTF.saturation_threshold", saturation_threshold);
		Prefs.set("PoissonNTF.bg_threshold", bg_threshold);
		Prefs.set("PoissonNTF.bg_choice", bg_choice);
		Prefs.set("PoissonNTF.subsamples", subsamples);
		Prefs.set("PoissonNTF.ndyes",ndyes);
		//initial spectra
		for(int dye=0; dye<ndyes; dye++)
		{			
			label="PoissonNTF.Dye_";
			label=label.concat(Integer.toString(dye+1));
			Prefs.set(label, spec_choice[dye]);
			label="PoissonNTF.DyeFixed_";
			label=label.concat(Integer.toString(dye+1));
			Prefs.set(label, spectra_fixed[dye]);
		}
	}

	/**
	 * Function that loops over all dyes and presents the user a save dialog
	 * to save each sprectrum
	 * @return
	 */
	public int saveEmissionSpectra(){
		//Dialog asking which spectra to save
		GenericDialog save_spectra=new GenericDialog("Save some spectra?");
		String [] labels=new String [ndyes];
		Boolean [] spectra_save=new Boolean [ndyes];
		for (int z = 0; z < ndyes; z++) {
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
			for (int z = 0; z < ndyes; z++) {
				spectra_save[z]=save_spectra.getNextBoolean();
			}
		}


		for (int dye = 0; dye < spectra_save.length; dye++) {
			if (spectra_save[dye])	//save requested spectra
			{
				SaveDialog sd=new SaveDialog("Save spectrum as ...", "spectrum"+new Integer(dye+1).toString(), ".emn");
				String dir=sd.getDirectory();
				String file=sd.getFileName();
				if (file==null) return -1;
				try{
					BufferedWriter out=new BufferedWriter(new FileWriter(dir+file));
					out.write("# Normalization: "+emnNorm+"\n");
					for (int emn = 0; emn < nemn; emn++){	//write wavelength in column 1, spectrum in column 2
						out.write(new Float(emission_wavelength.get(emn)).toString()+"\t"
								+new Float(S[emn][dye]/channel_width.get(emn)).toString()+"\t");
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

	
	/**
	 * Function that loops over all dyes and presents the user a save dialog
	 * to save each sprectrum
	 * @return
	 */
	public int saveExcitationSpectra(){
		//Dialog asking which spectra to save
		GenericDialog save_spectra=new GenericDialog("Save some spectra?");
		String [] labels=new String [ndyes];
		Boolean [] spectra_save=new Boolean [ndyes];
		for (int z = 0; z < ndyes; z++) {
			labels[z]= "Save excitation spectrum "+new Integer(z+1).toString()+"?";
			save_spectra.addCheckbox(labels[z], true);
		}
		save_spectra.showDialog();
		if (save_spectra.wasCanceled())
		{	
			return CANCELLED;
		}
		else
		{	//save result in boolean vector
			for (int z = 0; z < ndyes; z++) {
				spectra_save[z]=save_spectra.getNextBoolean();
			}
		}


		for (int dye = 0; dye < spectra_save.length; dye++) {
			if (spectra_save[dye])	//save requested spectra
			{
				SaveDialog sd=new SaveDialog("Save spectrum as ...", "spectrum"+new Integer(dye+1).toString(), ".exc");
				String dir=sd.getDirectory();
				String file=sd.getFileName();
				if (file==null) return -1;
				try{
					BufferedWriter out=new BufferedWriter(new FileWriter(dir+file));
					out.write("# Normalization: "+emnNorm+"\n");
					for (int exc = 0; exc < nexc; exc++){	//write wavelength in column 1, spectrum in column 2
						out.write(excitation_wavelength.get(exc)+"\t"
								+new Float(Q[exc][dye]).toString()+"\t");
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

	/**
	 * about dialog
	 */
	void showAbout() { 
		IJ.showMessage("About PoissonNTF...",
				"An ImageJ plugin for separating multiple dyes from an image stack. Version "+PoissonNTFversion
		); 
	} 

	private void calc_pinv()
	{
		int exc, dye1, dye2, emn, channel;
		double [][] STS=new double [ndyes][ndyes];
		double [][] STSch=new double [ndyes][ndyes];
		double [] y=new double [ndyes];

		for (dye1 = 0; dye1 < STS.length; dye1++) {
			for (dye2 = 0; dye2 < STS[dye1].length; dye2++) {
				STS[dye1][dye2]=0;
				for(exc=0; exc<nexc; exc++){
					for (emn = 0; emn < nchannels[exc]; emn++) {
					STS[dye1][dye2]+=Q[exc][dye1]*S[channeltoemission[exc][emn]][dye1]
					                 *Q[exc][dye2]*S[channeltoemission[exc][emn]][dye2];
					}
				}
			}
		}

		double sum,diff;
		int i,j,l, k1, k2;
		//cholesky decomposition
		for (i = 0; i < ndyes; i++) {
			sum = 0.0;
			for (j = 0; j < i; j++) {
				sum += STSch[i][j]*STSch[i][j];
			}
			diff = STS[i][i] - sum;
			STSch[i][i] = Math.sqrt(diff);
			for (l = i+1; l < ndyes; l++) {
				sum = 0.0;
				for (j = 0; j < i; j++) {
					sum += STSch[l][j]*STSch[i][j];
				}
				STSch[l][i] = (STS[l][i] - sum)/STSch[i][i];
				STSch[i][l]=0;
			}
		}
		//solving the lower diagonal system
		for (k1 = 0; k1 < ndyes; k1++) {
			for (k2 = 0; k2 < ndyes; k2++) y[k2]=0;
			y[k1]=1;
			for (i = 0; i < ndyes; i++) {
				sum = 0.0;
				for (j = 0; j < i; j++) {
					sum += STSch[i][j]*STS[j][k1];
				}
				STS[i][k1] = (y[i] - sum)/STSch[i][i];
			}	
		}
		//calculate the product of the inverse of the lower diagonal
		for (k1 = 0; k1 < ndyes; k1++) {
			for (k2 = k1; k2 <ndyes; k2++) {
				STSch[k1][k2] = 0.0;
				for (j = k2; j <ndyes; j++) {
					STSch[k1][k2] += STS[j][k1]*STS[j][k2];
				}
				STSch[k2][k1]=STSch[k1][k2];
			}	
		}
		//calculate the pseudo inverse
		for (dye1 = 0; dye1 < ndyes; dye1++) {
			channel=0;
			for(exc=0; exc<nexc; exc++){
				for (emn = 0; emn < nchannels[exc]; emn++) {
					pinvS[dye1][channel] = 0.0;
					for (dye2 = 0; dye2 < ndyes; dye2++) {
						pinvS[dye1][channel] += STSch[dye1][dye2]*S[channeltoemission[exc][emn]][dye2]*Q[exc][dye2];
					}
					channel++;
				}	
			}	
		}
	}

	
}
