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
 * OAS_ 
 * imagej plugin for OAS
 *
 * @author <a href="mailto:fabian.theis@helmholtz-muenchen.de">Fabian Theis</a>
 * @author <a href="mailto:neher@kitp.ucsb.edu">Richard Neher</a>
 */
public class OptimizeAcquisitionScheme_ implements PlugInFilter {
	String OASversion="0.1.1";

	//dimension of the images
	int ndyes=3;
	int nchannels=4;
	int nexc=2;
	int nemn;
	int nlaser=2;
	int totalchannels;
	int nemnpivots;
	int [][] channeltoemission;
	double []laserwavelength=null;
	double laserintensitystep=0.1;
	double []emnpivots=null;
	double gap=5;
	int [][]channels;
	int [][]optimalchannels;
	double [][]intensities;
	double [][]optimalintensities;
	double conditionnumber;
	double [][] noise;
	int lowerlambda=400;
	int upperlambda=700;
	int resolution=20;
	double bestconditionnumber_exc;
	double bestconditionnumber_emn;
	double[][] aobs_mask;
	int []first_channel; 
	//internal error codes
	int CANCELLED=-1;
	int NOROISELECTED=-2;
	int GOODINPUT=1;
	int BADVALUE=-3;
	double PMTNOISE=.05;
	//matrices that hold the objects and intermediates
	double[][] S;		//Emission spectra
	double[][] Strial;		//Emission spectra
	double[][] Q;		//Excitation spectra
	double[][] Qtrial;		//Excitation spectra
	double[][] pinvS;	//pseudo inverse of unfolded set of equations

	//auxillary quantities and parameters
	double [][][] Sexc;		//array holding the spectra from each excitation separately [nexc][nemm][ndyes]
	double E=1; 		//E=1 one-norm for spectra, E=2 two norm
	double intensitynothing=0.000001;
	boolean cancelled=false;
	boolean interrupted=false;
	boolean fixedchannels=false;
	boolean measure_blue=false;
	String [] spec_choice;			//array of strings with the choice of initial spectra
	PlotWindow plotEmission=null;			//window handle to display the spectra
	PlotWindow plotExcitation=null;			//window handle to display the spectra
	String config_file;
	double[] power; 	//buffer
	OptimizeAcquisitionProgress OASprogress;
	/**
	 * inital stuff
	 * @param arg
	 * @param imp
	 * @return
	 * @see ij.plugin.filter.PlugInFilter#setup(java.lang.String, ij.ImagePlus)
	 */
	public int setup(String arg, ImagePlus impin) {
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
			IJ.error("OAS requires ImageJ version 1.39t or greater!");
			return;
		}
		//user input of basic parameters
		if (parameterDialog()<0) return;
		/**
		 * normalize the segregation bias such that the contribution of the additional term to the cost
		 * function is of the order of that of a background pixel when segbias was 1
		 * the necessary magnitude of the segregation bias  depends on the overlap of the 
		 * spectra and the distribution of dye concentrations
		 */
		
		OASprogress=new OptimizeAcquisitionProgress(this);
		optimizeExcitations();
		optimizeChannels();
		IJ.showProgress(1,1);
		OASprogress.closeWindow();
		outputResults();
		//display results
		plotEmissionSpectra();
		plotExcitationSpectra();
	}

	/**
	 * Loop over all possible combinations of laser intensities for each excitation.
	 * the score of the setting is evalulated using the full resolution spectra matrix.
	 */
	void optimizeExcitations(){
		bestconditionnumber_exc=1e10;
		nemn=nemnpivots;
		//initialize intensities, calculate Q, S and the inverse
		initIntensities();
		initQtrial();
		calc_Sfull();
		calc_pinv();
		//use as best if the setting is non-singular
		if (!Double.isNaN(conditionnumber)) 
		{
			bestconditionnumber_exc=conditionnumber;
			saveexcitations();
		}
		//Loop over all possible
		//while loop terminates if next_excitation return false, indicating last setting
		while(next_excitationsetting()){
			//score new setting
			initQtrial();
			calc_Sfull();
			calc_pinv();
			if (conditionnumber<bestconditionnumber_exc  && !Double.isNaN(conditionnumber)){
				bestconditionnumber_exc=conditionnumber;
				saveexcitations();
			}						
		}
		
		//determine for each excitation the first channel that is on the red side of the 
		//reddest excitation
		if (measure_blue)
		{
			for (int exc=0; exc<nexc; exc++) first_channel[exc]=0;			
		}
		else{
			int exc, laser,emn;
			for (exc=0; exc<nexc; exc++)
			{
				double maxwl=0;
				for (laser=0; laser<nlaser; laser++)
				{
					if (optimalintensities[exc][laser]>intensitynothing && laserwavelength[laser]>maxwl) maxwl=laserwavelength[laser];
				}
				emn=0;
				do
				{
					emn++;
				}
				while (emnpivots[emn-1]<maxwl);
				first_channel[exc]=emn-1;
			}
		}
	}
	
	/**
	 * Use the previously determined best laserintensities and find the best
	 * channel setting for these lasers.
	 */
	void optimizeChannels()
	{
		bestconditionnumber_emn=1e10;
		nemn=nchannels;
		initQopt();
		initchannels();
		calc_pinv();
		if (!Double.isNaN(conditionnumber)) bestconditionnumber_emn=conditionnumber;
		savechannels();
		while(next_channelsetting() && cancelled==false){
			calc_Strial();
			calc_pinv();
			if (conditionnumber<bestconditionnumber_emn  && !Double.isNaN(conditionnumber)){
				bestconditionnumber_emn=conditionnumber;
				savechannels();
			}
	        IJ.showProgress(channels[0][nexc-1], nemnpivots-nchannels);
		}
	}

/**
 * set all channels in every excitation to 0 1 2 3 4 ....
 */
	void initchannels(){
		for(int exc=0; exc<nexc; exc++)
		{
			for(int ch=0; ch<nchannels+1; ch++)
			{
				channels[ch][exc]=first_channel[exc]+ch;
			}
		}
		calc_Strial();
	}
	
	/**
	 * set all intensities to laser0=1, all others 0 in every excitation
	 */
	void initIntensities(){
		for(int exc=0; exc<nexc; exc++)
		{
			intensities[exc][0]=1;
			for(int laser=1; laser<nlaser; laser++)
			{
				intensities[exc][laser]=0;
			}
		}
		calc_Strial();
	}
	
	/**
	 * save the current channel setting as the best-up-to now
	 */
	void savechannels(){
		for(int exc=0; exc<nexc; exc++)
		{
			for(int ch=0; ch<nchannels+1; ch++)
			{
				optimalchannels[ch][exc]=channels[ch][exc];
			}
		}
	}

	/**
	 * save the current setting laser intensities as the best-up-to now
	 */
	void saveexcitations(){
		for(int exc=0; exc<nexc; exc++)
		{
			for(int laser=0; laser<nlaser; laser++)
			{
				optimalintensities[exc][laser]=intensities[exc][laser];
			}
		}
		
	}
	
	/**
	 * function that increments the channel setting by starting at the last channel and 
	 * moving it through the available range, if that is exhausted, it increments the next
	 * to last while resetting the channels to the right, and so on for channels to left. 
	 * if all channels are in their right most position, it moves on the next excitation 
	 * resets the previous excitations
	 * @return true if valid setting is found, false, if the very last setting is encountered
	 */
	boolean next_channelsetting()
	{
		int exc=0, exc2;
		boolean check=false;
		int ch=nchannels;	//start with last channel
		do
		{
			if (channels[ch][exc]>(nemnpivots-(nchannels-ch)-2)) ch--;	//channel cannot be moved, try next channel
			else {channels[ch][exc]++; check=true;}	//channel can be moved, exit loop
			if (fixedchannels==false && ch<0 && exc<(nexc-1)) {ch=nchannels; exc++;}	//if all channels have been tried, move on the next excitation
		}while(check==false && (ch>=0));
		//if a channel has been moved, reset the channel to the right and all prior excitations
		if (ch>=0) 
		{
			for (++ch;ch<nchannels+1;ch++) channels[ch][exc]=channels[ch-1][exc]+1;
			for (exc2=0; exc2<exc; exc2++)
			{
				for (ch=0;ch<nchannels+1;ch++) channels[ch][exc2]=first_channel[exc2]+ch;
			}
		}
		
		//if the same channels are to be used in each excitation, copy
		if (fixedchannels)
		{
			for (exc=1; exc<nexc; exc++){
				for (ch=0; ch<nchannels+1; ch++) channels[ch][exc]=channels[ch][0];
			}
		}
		return check;
	}
	
	
	/**
	 * loop over laser intensities, such that the intensities sum to one at each excitation.
	 * @return
	 */
	boolean next_excitationsetting()
	{
		boolean check=false;
		int exc=0, laser=0;
		do{
			check=incrementExcitation(exc, 1, 0);
			if (check==false) 
			{
				intensities[exc][0]=1;
				for (laser=1; laser<nlaser; laser++) intensities[exc][laser]=0;
				exc++; 
				if (exc<nexc) check=incrementExcitation(exc, 1, 0);
			}
		}while(exc<nexc && check==false);
		return check;
	}
	
	/**
	 * recursive function that loops over all possible ways to distribute remainder over 
	 * lasers laser..nlaser-1
	 * @param exc	exctitation
	 * @param remainder	the remaingn excitiation allowance
	 * @param laser	specifies the lasers among which the remainder is shared
	 * @return
	 */
	boolean incrementExcitation(int exc, double remainder, int laser)
	{
		boolean check;
		if (laser>=nlaser-1 || remainder<-0.01*laserintensitystep) {remainder=0; return false;}
		else  check=incrementExcitation(exc, remainder-intensities[exc][laser], laser+1);
		
		if (check) return true;
		else if (intensities[exc][laser]-laserintensitystep>-0.01*laserintensitystep) {
			intensities[exc][laser]-=laserintensitystep;
			intensities[exc][laser+1]=remainder-intensities[exc][laser];
			for (int l=laser+2; l<nlaser; l++) intensities[exc][l]=0;
			return true;
		}
		else return false;
	}
	
	/**
	 * calculate the mixing matrix based on the current channels and the 
	 * optimal exciation settings
	 */
	void calc_Strial()
	{
		int exc, ch, emn, dye;
		for(exc=0; exc<nexc; exc++)
		{
			for(ch=0; ch<nchannels; ch++)
			{
				noise[exc][ch]=(channels[ch+1][exc]-channels[ch][exc])*PMTNOISE;	//some extra variance for each channel
				for(dye=0; dye<ndyes; dye++)
				{
					Sexc[exc][ch][dye]=0;
					for (emn=channels[ch][exc]; emn<channels[ch+1][exc]; emn++)
					{
						Sexc[exc][ch][dye]+=S[emn][dye];
					}
					Sexc[exc][ch][dye]*=Qtrial[exc][dye];
					noise[exc][ch]+=Sexc[exc][ch][dye];	//the variance is proportional to the signal
				}
				noise[exc][ch]=Math.sqrt(noise[exc][ch]);
				for(dye=0; dye<ndyes; dye++)
				{
					Sexc[exc][ch][dye]/=noise[exc][ch];	//rescale the matrix such that each column has unit variance
				}
			}
		}
	}
	
	/**
	 * calculate the mixing matrix for the full spectral resolution case
	 */
	void calc_Sfull()
	{
		int exc, ch, dye,laser;
		for(exc=0; exc<nexc; exc++)
		{
			for(ch=0; ch<nemnpivots; ch++)
			{
				noise[exc][ch]=PMTNOISE;
				for(dye=0; dye<ndyes; dye++)
				{
					Sexc[exc][ch][dye]=S[ch][dye]*Qtrial[exc][dye];
					noise[exc][ch]+=Sexc[exc][ch][dye];
					for (laser=0; laser<nlaser; laser++)
					{
						if (intensities[exc][laser]>0) Sexc[exc][ch][dye]*=aobs_mask[laser][ch];
					}
				}
				noise[exc][ch]=Math.sqrt(noise[exc][ch]);
				for(dye=0; dye<ndyes; dye++)
				{
					Sexc[exc][ch][dye]/=noise[exc][ch];	//rescale the matrix such that each channel has unit variance
				}
			}
		}
	}

	/**
	 * initialize the matrix Qtrial with current intensities
	 * to be used in calc_Sfull
	 */
	void initQtrial()
	{
		for (int exc=0; exc<nexc; exc++)
		{
			for (int dye=0; dye<ndyes; dye++){
				Qtrial[exc][dye]=0;
				for (int laser=0; laser<nlaser; laser++) Qtrial[exc][dye]+=intensities[exc][laser]*Q[laser][dye];
			}
		}
	}

	/**
	 * initialize the matrix Qtrial with the optimal intensities
	 */
	void initQopt()
	{
		for (int exc=0; exc<nexc; exc++)
		{
			for (int dye=0; dye<ndyes; dye++){
				Qtrial[exc][dye]=0;
				for (int laser=0; laser<nlaser; laser++) Qtrial[exc][dye]+=optimalintensities[exc][laser]*Q[laser][dye];
			}
		}
		
	}
	
	/**
	 * Print optimal solution to the imageJ results window
	 */
	void outputResults(){
		String out;
		IJ.write("Score of the optimal acquisition setting with full spectral resolution: ");
		IJ.write(Double.toString(bestconditionnumber_exc));
		IJ.write("Score of the optimal acquisition setting: ");
		IJ.write(Double.toString(bestconditionnumber_emn));
		for (int exc=0; exc<nexc; exc++){
			IJ.write("\nExcitation "+Integer.toString(exc+1)+":");
			for (int laser=0; laser<nlaser; laser++) {
				out="Laser ";
				out=out.concat(IJ.d2s(laserwavelength[laser],2)+"nm: "+IJ.d2s(optimalintensities[exc][laser],2));
				IJ.write(out);
			}
			for (int ch=0; ch<nchannels; ch++){
				IJ.write("Channel "+Integer.toString(ch+1)+": "+
						Double.toString(emnpivots[optimalchannels[ch][exc]])+"-"+
						Double.toString(emnpivots[optimalchannels[ch+1][exc]])+"nm");
			}			
			IJ.write("\n");
		}
	}
	
	/**
	 * Prompt user for input, get background and start spectra
	 */
	private int parameterDialog()
	{
		GenericDialog nsources= new GenericDialog("OAS");
		nsources.addNumericField("Number of Sources", ndyes, 0);
		nsources.addNumericField("Number of Lasers", nlaser, 0);
		nsources.showDialog();
		if (nsources.wasCanceled())
		{	
			return CANCELLED;
		}
		else
		{
			ndyes=(int)nsources.getNextNumber();
			if (ndyes<2 || ndyes>10) 
			{
				IJ.error("Bad number of Sources!");
				return -1;
			}
			nlaser=(int)nsources.getNextNumber();
			if (nlaser<1 || ndyes>10) 
			{
				IJ.error("Bad number of lasers!");
				return -1;
			}

			allocate_ndyes_dependent();
			//
			//get preferences from IJ preferences file
			getPrefs();

			GenericDialog OAS_dialog=new GenericDialog("Optimize Acquisition Scheme");
			OAS_dialog.addNumericField("Lambda min [nm]: ", lowerlambda, 0);
			OAS_dialog.addNumericField("Lambda max [nm]: ", upperlambda, 0);
			OAS_dialog.addNumericField("Lambda steps [nm]: ", resolution, 0);
			OAS_dialog.addNumericField("Number of exc: ", nexc, 0);
			OAS_dialog.addNumericField("Number of channels: ", nchannels, 0);
			OAS_dialog.addNumericField("Laser intensity steps [%]: ", laserintensitystep*100, 0);
			OAS_dialog.addCheckbox("Same channels in all excitations?", fixedchannels);
			OAS_dialog.addCheckbox("Allow excitations within emission?", measure_blue);
			add_lasers(OAS_dialog);
			add_spectra_choice(OAS_dialog);
			OAS_dialog.setOKLabel("Optimize!");

			OAS_dialog.showDialog();
			if (OAS_dialog.wasCanceled())
			{	
				return CANCELLED;
			}
			else
			{
				//RUNTIME
				//number of iterations
				lowerlambda=(int) OAS_dialog.getNextNumber();
				upperlambda=(int) OAS_dialog.getNextNumber();
				resolution=(int) OAS_dialog.getNextNumber();
				nexc=(int) OAS_dialog.getNextNumber();
				nemnpivots=(upperlambda-lowerlambda)/resolution;
				nchannels=(int) OAS_dialog.getNextNumber();
				laserintensitystep=OAS_dialog.getNextNumber()/100;
				fixedchannels=OAS_dialog.getNextBoolean();
				measure_blue=OAS_dialog.getNextBoolean();
				allocate_nexc_dependent();			
				allocate_nemn_dependent();
				read_laser_wavelength(OAS_dialog);
				read_spectra_choice(OAS_dialog);
			}
		}
		setPrefs();
		return 0;
	}

	/**
	 * init S with Gaussians
	 */
	private void gauss_spectra() {
		int emn;
		int dye;
		//spectra - Gauss shape
		for (dye=0;dye<ndyes;dye++) {		//loop over dyes
			float sigma = (nemnpivots-1f)/(ndyes+1f)+.1f;	//variance of gauss spectra
			for (emn=0;emn<nemnpivots;emn++) 				//loop over channels and assign spectra
			{
				S[emn][dye]= (1f/(sigma*Math.sqrt(2*Math.PI))*
						Math.exp(-.5*Math.pow(emn-1-((dye+1)*nemnpivots)/(ndyes+1f),2)
								/Math.pow(sigma,2))*(resolution));
			}
		}
	}


	
	/**
	 * add the available laser wavelength to the dialog
	 */
	private int add_lasers(GenericDialog choose_lasers)
	{
		String label;
		for(int laser=0; laser<nlaser; laser++){
			label="Laserwavelength ";
			label=label.concat(Integer.toString(laser+1));
			label=label.concat(" [nm]");			
			choose_lasers.addNumericField(label, laserwavelength[laser], 0);
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
		FilenameFilter specemn_file_filter = new FilenameFilter() {
			public boolean accept(File dir, String name) {
				return name.endsWith(".emn");
			}
		};
		//filter files that end on .exc (for exc spectra)
		FilenameFilter specexc_file_filter = new FilenameFilter() {
			public boolean accept(File dir, String name) {
				return name.endsWith(".exc");
			}
		};
		//produce a list of files that end on emn in directory /plugins/SpectraLibrary
		String curDir = System.getProperty("user.dir");
		curDir=curDir.concat("/plugins/SpectraLibrary");
		File spec_directory=new File(curDir);
		String[] spectra_list=spec_directory.list(specemn_file_filter);
		String[] items;
		String default_choice;
		if (spectra_list!=null)
		{
			items=new String [spectra_list.length + 2];
		}else{
			items=new String [2];			
		}
		//the choices for the spectra + the content of the directory /plugins/SpectraLibrary
		items[0]="manually";
		items[1]="--------";

		for (int i = 2; i < items.length; i++) {
			items[i]=spectra_list[i-2];
		}

		String label;
		for (int dye = 0; dye < ndyes; dye++) {
			default_choice=items[0];
			if (spec_choice[dye]!=null)
			{
				for (int i=0; i<items.length; i++)
				{
					if (spec_choice[dye].equals(items[i])) 	default_choice=items[i];
				}
			}
			label="Emission spectrum dye ";
			label=label.concat(Integer.toString(dye+1));
			choose_spectra.addChoice(label, items, default_choice);
		}		
		
		//produce a list of files that end on exc in directory /plugins/SpectraLibrary
		spectra_list=spec_directory.list(specexc_file_filter);
		if (spectra_list!=null)
		{
			items=new String [spectra_list.length + 2];
		}else{
			items=new String [2];			
		}
		//the choices for the spectra + the content of the directory /plugins/SpectraLibrary
		items[0]="manually";
		items[1]="--------";

		for (int i = 2; i < items.length; i++) {
			items[i]=spectra_list[i-2];
		}
		for (int dye = 0; dye < ndyes; dye++) {
			default_choice=items[0];
			for (int i=0; i<items.length; i++)
			{
				if (spec_choice[dye+ndyes]!=null)
					if (spec_choice[dye+ndyes].equals(items[i])) 	default_choice=items[i];
			}
			label="Exc. eff. dye ";
			label=label.concat(Integer.toString(dye+1));
			choose_spectra.addChoice(label, items, default_choice);
		}
		return 0;

	}

	/**
	 * handle the choice of spectra from the parameter dialog
	 *
	 */
	private int read_spectra_choice(GenericDialog choose_spectra)
	{
		//the three choices for the background
		String[] items=new String [2];
		String curDir = System.getProperty("user.dir");
		curDir=curDir.concat("/plugins/SpectraLibrary");
		
		gauss_spectra();
		
		//the choices for the spectra + the content of the directory /plugins/SpectraLibrary
		items[0]="manually";
		items[1]="--------";
		//deal with background first, it is needed for the determination of the spectra later

		//get choice for each spectrum and assign spectra
		for (int dye = 0; dye < ndyes; dye++) {
			spec_choice[dye]=choose_spectra.getNextChoice();
			if (spec_choice[dye].equals(items[0])) enterSpectrum(dye);		//Manually
			else if (spec_choice[dye].equals(items[1])) {}
			else readSpectrum(dye, curDir, spec_choice[dye]);				//from file
		}

		//get choice for each spectrum and assign spectra
		for (int dye = 0; dye < ndyes; dye++) {
			spec_choice[dye+ndyes]=choose_spectra.getNextChoice();
			if (spec_choice[dye+ndyes].equals(items[0])) enterExcitationSpectrum(dye);		//Manually
			else if (spec_choice[dye+ndyes].equals(items[1])) {}
			else readExcitationSpectrum(dye, curDir, spec_choice[dye+ndyes]);				//from file
		}

		
		//Normalize the spectra
		for (int dye=0;dye<ndyes;dye++) {				//for each dye, calculate the 'E'-norm of the spectrum
			power[dye]=0.0;
			for (int emn=0;emn<nemnpivots;emn++)
				power[dye]+=Math.pow(S[emn][dye],E);
			power[dye]=(float) Math.pow(power[dye], 1f/E);
		}
		for (int dye=0;dye<ndyes;dye++) {			//rescale S
			for (int emn=0;emn<nemnpivots;emn++)
			{
				S[emn][dye]/=power[dye];
			}
		}
		return 0;
	}
	
	/**
	 * handle the choice of spectra from the parameter dialog
	 *
	 */
	private int read_laser_wavelength(GenericDialog choose_spectra)
	{
		//the three choices for the background
		for (int laser = 0; laser < nlaser; laser++) {
			laserwavelength[laser]=choose_spectra.getNextNumber();
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
			int emn=0;
			double step, result;	//integration step
			for (emn = 0; emn < nemnpivots; emn++) 
			{
				step=resolution/20.0;	//integration step
				result=0;
				for (double lambdax=emnpivots[emn]-0.5*resolution;lambdax<emnpivots[emn]+0.5*resolution; lambdax+=step)
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
	 * Read spectrum form file. 
	 * @param dye
	 * @param dir
	 * @param spec
	 */
	private void readExcitationSpectrum(int dye, String dir, String spec)
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
			int exc=0;
			for (exc = 0; exc < nlaser; exc++) 
			{
				Q[exc][dye]=interpolate(lambda, emission, laserwavelength[exc]);
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
		if (lambda_0<=lambda.get(0)) //value outside the interval -> extrapolate
		{
			double lambda1=lambda.get(0), lambda2=lambda.get(1), F1=F.get(0), F2=F.get(1);
			result=F1+(F2-F1)/(lambda2-lambda1)*(lambda_0-lambda1);
			result=Math.min(F1, F2);
			result=0;
		}
		else if (lambda_0>=lambda.get(lambda.size()-1))  //value outside the interval -> extrapolate
		{
			double lambda1=lambda.get(lambda.size()-2), lambda2=lambda.get(lambda.size()-1), F1=F.get(lambda.size()-2), F2=F.get(lambda.size()-1);
			result=F1+(F2-F1)/(lambda2-lambda1)*(lambda_0-lambda1);
			result=0;
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
		for (int emn=0; emn<nemnpivots; emn++)
		{	//add numeric field with the current spectrum as default value
			enter_spectra.addNumericField("Wavelength "+Double.toString(emnpivots[emn])+"nm", S[emn][dye], 1);
		}
		enter_spectra.showDialog();
		if (enter_spectra.wasCanceled())
		{	
			return CANCELLED;
		}
		else
		{
			for (int emn=0; emn<nemnpivots; emn++)
			{	//get user input
				S[emn][dye]=enter_spectra.getNextNumber();
			}
		}
		return GOODINPUT;
	}

	/**
	 * Open a dialog with a numeric field for each channel such that the user can enter 
	 * the spectrum of the dye 
	 * @param dye
	 * @return
	 */
	private int enterExcitationSpectrum(int dye)
	{
		String msg;
		msg="Exc. eff. of dye ";
		msg=msg.concat(Integer.toString(dye+1));
		GenericDialog enter_spectra=new GenericDialog(msg);
		for (int exc=0; exc<nlaser; exc++)
		{	//add numeric field with the current spectrum as default value
			enter_spectra.addNumericField("Laser: "+Double.toString(laserwavelength[exc])+"nm", Q[exc][dye], 1);
		}
		enter_spectra.showDialog();
		if (enter_spectra.wasCanceled())
		{	
			return CANCELLED;
		}
		else
		{
			for (int exc=0; exc<nlaser; exc++)
			{	//get user input
				Q[exc][dye]=enter_spectra.getNextNumber();
			}
		}
		return GOODINPUT;
	}



	/**
	 * generate new plot for the spectra
	 * @return 
	 * 
	 */
	public void plotEmissionSpectra() {
		int dye, emn;
		double [][] spec=new double [nemnpivots][ndyes];
		for (dye=0;dye<ndyes;dye++) {
			for (emn=0;emn<nemnpivots;emn++) {
				spec[emn][dye]=S[emn][dye];
			}
		}
		//colors for the different spectra
		Color[] colors = {Color.blue,Color.green, Color.red, Color.cyan, Color.gray, Color.darkGray};
		Plot plot = new Plot("OAS emission spectra","wave length [nm]","intensity",(float[])null,(float[])null, PlotWindow.LINE);

		if (plotEmission==null){
			plotEmission=new PlotWindow("OAS emission spectra","wave length [nm]","intensity",(float[])null,(float[])null);
		}
		float[] rowS = new float[nemnpivots];
		float[] lambdas=new float [nemnpivots];
		for (emn = 0; emn < lambdas.length; emn++) {
			lambdas[emn]=(float) (emnpivots[emn]); //channel wavelength as x-axis
		}
		plot.setLimits(emnpivots[0], emnpivots[nemnpivots-1], 0, 1);
		int c=0; float maxS=0;
		//determine maximum of all spectra
		for (dye=0;dye<ndyes;dye++) {
			for (emn=0;emn<nemnpivots;emn++) {
				if (spec[emn][dye]/resolution>maxS) 
					maxS=(float) ((float) spec[emn][dye]/resolution);
			}
		}
		//loop over spectra and plot each
		for (dye=0;dye<ndyes;dye++) {
			for (emn=0;emn<nemnpivots;emn++) {
				rowS[emn]=(float) (spec[emn][dye]/resolution/maxS);
			}
			plot.setLineWidth(2);            
			plot.setColor(colors[c]);
			plot.addPoints(lambdas, rowS, PlotWindow.LINE);
			//plot initial spectra only if not kept fixed
			if (c++>colors.length) c=0;
		}
		double [] channelboundary_x=new double [2];
		double [] channelboundary_y=new double [2];		
		for (int exc=0; exc<nexc; exc++)
		{
			channelboundary_y[0]=1.0*exc/nexc;
			channelboundary_y[1]=1.0*(exc+1)/nexc;
			plot.setColor(colors[exc]);
			for (int ch=0; ch<nchannels+1; ch++)
			{
				channelboundary_x[0]=emnpivots[optimalchannels[ch][exc]];
				channelboundary_x[1]=emnpivots[optimalchannels[ch][exc]];
				plot.setLineWidth(2);            
				plot.addPoints(channelboundary_x, channelboundary_y, PlotWindow.LINE);
			}
		}
		plot.draw();
		plotEmission.drawPlot(plot);
	}

	/**
	 * generate new plot for the spectra
	 * @return 
	 * 
	 */
	public void plotExcitationSpectra() {
		int dye, exc;
		double [][] spec=new double [nlaser][ndyes];
		for (dye=0;dye<ndyes;dye++) {
			for (exc=0;exc<nlaser;exc++) {
				spec[exc][dye]=Q[exc][dye];
			}
		}
		//colors for the different spectra
		Color[] colors = {Color.blue,Color.green, Color.red, Color.cyan, Color.gray, Color.darkGray};
		Plot plot = new Plot("PoissonNMF excitation spectra","Excitation","intensity",(float[])null,(float[])null, PlotWindow.LINE);

		if (plotExcitation==null){
			plotExcitation=new PlotWindow("OAS excitation spectra","Excitation","intensity",(float[])null,(float[])null);
		}
		float[] rowS = new float[nlaser];
		float[] lambdas=new float [nlaser];
		for (exc = 0; exc < lambdas.length; exc++) {
			lambdas[exc]=(float) laserwavelength[exc];//channel wavelength as x-axis
		}
		int c=0; float maxS=0;
		//determine maximum of all spectra
		for (dye=0;dye<ndyes;dye++) {
			for (exc=0;exc<nlaser;exc++) {
				if (spec[exc][dye]>maxS) 
					maxS=(float) spec[exc][dye];
			}
		}
		plot.setLimits(lambdas[0]-0.5, lambdas[nlaser-1]+0.5, 0, maxS);
		//loop over spectra and plot each
		for (dye=0;dye<ndyes;dye++) {
			for (exc=0;exc<nlaser;exc++) {
				rowS[exc]=(float) spec[exc][dye];
			}
			plot.setLineWidth(2);            
			plot.setColor(colors[c]);
			plot.addPoints(lambdas, rowS, PlotWindow.LINE);
			if (c++>colors.length) c=0;
		}
		plot.draw();
		plotExcitation.drawPlot(plot);
	}


	private void allocate_nexc_dependent()
	{
		//Vectors that save the emission channel boundaries for each excitation.
		Sexc = new double [nexc][nemnpivots][ndyes];
		noise = new double [nexc][nemnpivots];
		Q = new double [nlaser][ndyes];
		Qtrial = new double [nexc][ndyes];
		intensities = new double [nexc][nlaser];
		optimalintensities = new double [nexc][nlaser];
		for(int exc=0; exc<nexc; exc++) totalchannels+=nemnpivots;
		pinvS=new double [ndyes][totalchannels];
		first_channel=new int [nexc];
		for (int exc=0; exc<nexc; exc++) first_channel[exc]=0;
	}

	private void allocate_ndyes_dependent() {
		spec_choice=new String [2*ndyes];
		laserwavelength=new double [nlaser];
	}


	private void allocate_nemn_dependent()
	{
		S = new double [nemnpivots][ndyes];
		power=new double [ndyes];
		channeltoemission= new int [nexc][nchannels];
		optimalchannels=new int [nchannels+1][nexc];
		channels=new int [nchannels+1][nexc];
		emnpivots= new double [nemnpivots];
		for(int emn=0; emn<nemnpivots; emn++) emnpivots[emn]=lowerlambda+emn*resolution;
		aobs_mask= new double [nlaser][nemnpivots];
		for (int laser=0; laser<nlaser; laser++)
		{
			for(int ch=0; ch<nemnpivots; ch++)
			{
				if (emnpivots[ch]<laserwavelength[laser] && laserwavelength[laser]<emnpivots[ch]+resolution)
				{
					aobs_mask[laser][ch]=(resolution-gap)/resolution;
				}
				else aobs_mask[laser][ch]=1.0;
			}			
		}

	}

	/**
	 * Function that reads the previous settings from the imageJ preferences file
	 **/
	private void getPrefs() {
		String label;
		nexc=(int) Prefs.get("OAS.nexc",nexc);
		nchannels=(int) Prefs.get("OAS.nchannels",nchannels);
		laserintensitystep=Prefs.get("OAS.intensitysteps", laserintensitystep);
		resolution=(int) Prefs.get("OAS.resolution", resolution);
		lowerlambda=(int) Prefs.get("OAS.lambdamin", lowerlambda);
		upperlambda=(int) Prefs.get("OAS.lambdamax", upperlambda);
		//choice of initial spectra and decision to keep some spectra fixed
		if ((int)Prefs.get("OAS.ndyes",0)==ndyes)
		{
			for(int dye=0; dye<2*ndyes; dye++)
			{			
				label="OAS.Dye_";
				label=label.concat(Integer.toString(dye+1));
				spec_choice[dye]=Prefs.get(label, "none");
			}
		}
		for(int laser=0; laser<nlaser; laser++)
		{
			label="OAS.Laser_";
			label=label.concat(Integer.toString(laser+1));
			laserwavelength[laser]=Prefs.get(label, laserwavelength[laser]);
		}
	}

	/**
	 * Save the parameters to the imageJ preferences file: one-to-one correspondence
	 * to getPrefs() above
	 */
	private void setPrefs() {
		String label;
		Prefs.set("OAS.ndyes",ndyes);
		Prefs.set("OAS.nchannels",nchannels);
		Prefs.set("OAS.nexc",nexc);
		Prefs.set("OAS.intensitysteps", laserintensitystep);
		Prefs.set("OAS.resolution", resolution);
		Prefs.set("OAS.lambdamin", lowerlambda);
		Prefs.set("OAS.lambdamax", upperlambda);
		//initial spectra
		for(int dye=0; dye<2*ndyes; dye++)
		{			
			label="OAS.Dye_";
			label=label.concat(Integer.toString(dye+1));
			Prefs.set(label, spec_choice[dye]);
		}
		for(int laser=0; laser<nlaser; laser++)
		{
			label="OAS.Laser_";
			label=label.concat(Integer.toString(laser+1));
			Prefs.set(label, laserwavelength[laser]);
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
					for (int emn = 0; emn < nchannels; emn++){	//write wavelength in column 1, spectrum in column 2
						out.write(new Float(emnpivots[emn]).toString()+"\t"
								+new Float(S[emn][dye]/resolution).toString()+"\t");
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
					for (int exc = 0; exc < nexc; exc++){	//write wavelength in column 1, spectrum in column 2
						out.write(new Float(emnpivots[exc]).toString()+"\t"
								+new Float(Q[exc][dye]/resolution).toString()+"\t");
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
		IJ.showMessage("About OAS...",
				"An ImageJ plugin for separating multiple dyes from an image stack. Version "+OASversion
		); 
	} 

	private void calc_pinv()
	{
		int exc, dye1, dye2, ch;
		double [][] STS=new double [ndyes][ndyes];
		double [][] STSch=new double [ndyes][ndyes];
		//double [] y=new double [ndyes];

		for (dye1 = 0; dye1 < STS.length; dye1++) {
			for (dye2 = 0; dye2 < STS[dye1].length; dye2++) {
				STS[dye1][dye2]=0;
				for(exc=0; exc<nexc; exc++){
					for (ch = 0; ch < nemn; ch++) {
					STS[dye1][dye2]+=Sexc[exc][ch][dye1]*Sexc[exc][ch][dye2];
					}
				}
			}
		}

		double sum,diff;
		int i,j,l, k1; // k2;
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
		conditionnumber=-1;
		for (k1 = 0; k1 < ndyes; k1++) conditionnumber*=STSch[k1][k1];
		/*
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
		}*/
/*		conditionnumber=0;
		for (k1 = 0; k1 < ndyes; k1++) conditionnumber+=STSch[k1][k1];*/
			
/*		//calculate the pseudo inverse
		for (dye1 = 0; dye1 < ndyes; dye1++) {
			channel=0;
			for(exc=0; exc<nexc; exc++){
				for (emn = 0; emn < nemn; emn++) {
					pinvS[dye1][channel] = 0.0;
					for (dye2 = 0; dye2 < ndyes; dye2++) {
						pinvS[dye1][channel] += STSch[dye1][dye2]*Sexc[exc][emn][dye2];
					}
					channel++;
				}	
			}	
		}
*/
	}

	
}
