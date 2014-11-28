import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;

import ij.*;
import ij.gui.*;

	public class PoissonNMFspectracorrection implements MouseListener{
	int ndyes=3;
	double [][] G;
	double [][] Ginv;
	double [] base1;
	double [] base2;
	double [] base3;
	double alpha;
	double xscale, xoffset;
	double yscale, yoffset;
	int [] dyes;
	ImagePlus triangle_image;
	ImageCanvas canvas;
	int maxpix=400;
	float [][] conc_hist;
	int [][] spec_triangle;
	int [][] positivity_boundary;
	PoissonNMF_ parent;

	int spectraChosen;
	boolean pressed;
	
	public PoissonNMFspectracorrection(PoissonNMF_ pnmf, int [] dyes_in)
	{
		parent=pnmf;	
		if (dyes_in.length!=3) return;
		dyes=new int [dyes_in.length];
		for (int i = 0; i < dyes_in.length; i++) {
			dyes[i] = dyes_in[i];
		}
		//		xmax=((double)index_x(1,0))/maxpix*1.2;	ymax=((double)index_y(0,1))/maxpix*1.2;
		//		xmin=-((double)index_x(1,0))/maxpix*0.1;	ymin=-((double)index_y(0,1))/maxpix*0.1;
		//double [] pixconc=new double [pnmf.r];
		conc_hist=new float [maxpix][maxpix];
		spec_triangle=new int [maxpix][maxpix];
		positivity_boundary=new int [maxpix][maxpix];
		triangle_image=NewImage.createByteImage("Data simplex", maxpix, maxpix, 1, 0);

		default_values();
		calc_and_draw();
	}

	private void calc_and_draw(){
		calc_base();
		calc_positivity_region();
		calc_spec_triangle();
		calc_conc_histogram();
		updateImage();		
	}
	
	
	private void calc_positivity_region() {
		boolean N=true, S=true, E=true, W=true;
		double muN=0, muS=1, muE=0, muW=1;
		double mu1, mu2;
		double [] spec= new double [parent.n];
		int test1, test2;
		//find the a rectangle such that the entire boundary is in the negative spectra region
		while (N || S || E || W)
		{
			//do the same for all four boundaries
			if (N==true) 
			{	
				muN+=0.01;	//push boundary out
				mu2=-muN; 
				test1=0;
				for (mu1=-muE; mu1<=muW; mu1+=0.01)	//walk along the boundary
				{	
					test2=0;
					for(int lambda=0; lambda<spec.length; lambda++)
					{	//increment test2 if any channel is negative
						if (base3[lambda]+mu1*base1[lambda]+mu2*base2[lambda]<0) test2++;
					}
					if (test2==0) test1++;	//increment test1 if all channels were positive
				}
				if (test1==0)	//if test1 didn't get incremented, all pixels on the boundary violated positivity 
					N=false;
			}
			if (S==true) 
			{
				muS+=0.01;
				mu2=muS; 
				test1=0;
				for (mu1=-muE; mu1<=muW; mu1+=0.01)
				{	
					test2=0;
					for(int lambda=0; lambda<spec.length; lambda++)
					{
						if (base3[lambda]+mu1*base1[lambda]+mu2*base2[lambda]<0) test2++;
					}
					if (test2==0) test1++;
				}
				if (test1==0)
					S=false;
			}
			if (E==true) 
			{
				muE+=0.01;
				mu1=-muE; 
				test1=0;
				for (mu2=-muN; mu2<=muS; mu2+=0.01)
				{	
					test2=0;
					for(int lambda=0; lambda<spec.length; lambda++)
					{
						if (base3[lambda]+mu1*base1[lambda]+mu2*base2[lambda]<0) test2++;
					}
					if (test2==0) test1++;
				}
				if (test1==0)
					E=false;
			}

			if (W==true) 
			{
				muW+=0.01;
				mu1=muW; 
				test1=0;
				for (mu2=-muN; mu2<=muS; mu2+=0.01)
				{	
					test2=0;
					for(int lambda=0; lambda<spec.length; lambda++)
					{
						if (base3[lambda]+mu1*base1[lambda]+mu2*base2[lambda]<0) test2++;
					}
					if (test2==0) test1++;
				}
				if (test1==0)
					W=false;
			}
		}
		//calculate the the offset and scale parameter used when calculating indices such 
		//that the entire positive domain fits into the image
		yscale=G[1][1]*Math.sin(alpha)*(muS+muN);
		yoffset=-G[1][1]*Math.sin(alpha)/yscale*muN;
		xscale=(G[0][0]*(muW+muE)+G[1][1]*Math.cos(alpha)*(muS+muN));
		xoffset=-(G[0][0]*muE+G[1][1]*Math.cos(alpha)*muN)/xscale;
	
		//make a mask that masks out the negative forbidden regions
		for (int i = 0; i < positivity_boundary.length; i++) {
			for (int j = 0; j < positivity_boundary.length; j++) {
				positivity_boundary[i][j]=0;
			}
		}
		for (int i = 0; i < positivity_boundary.length; i++) {
			for (int j = 0; j < positivity_boundary.length; j++) {
				mu1=invindex_x(i,j);
				mu2=invindex_y(i,j);
				test2=0;
				for(int lambda=0; lambda<spec.length; lambda++)
				{
					if (base3[lambda]+mu1*base1[lambda]+mu2*base2[lambda]<0) test2++;
				}
				if (test2>0) positivity_boundary[i][j]=50;

			}
		}		

	}

	private void calc_spec_triangle() {
		int i1, i2;
		for (int i = 0; i < spec_triangle.length; i++) {
			for (int j = 0; j < spec_triangle.length; j++) {
				spec_triangle[i][j]=0;
			}
		}		
		for (double x=0; x<=1; x+=0.001)
		{
			i1=index_x(x, 0);
			i2=index_y(x, 0);
			if (i2>=0 && i2<maxpix && i1>=0 && i1<maxpix) spec_triangle[i1][i2]=100;
			i1=index_x(0, x);
			i2=index_y(0, x);
			if (i2>=0 && i2<maxpix && i1>=0 && i1<maxpix) spec_triangle[i1][i2]=100;
			i1=index_x(1-x, x);
			i2=index_y(1-x, x);
			if (i2>=0 && i2<maxpix && i1>=0 && i1<maxpix) spec_triangle[i1][i2]=100;
		}
	}

	private void calc_conc_histogram()
	{
		int i1, i2;
		for (i1=0; i1<conc_hist.length; i1++){
			for (i2=0; i2<conc_hist.length; i2++){
				conc_hist[i1][i2]=0;
			}
		}
		double pixnorm, scalarprod1, scalarprod2, s1, s2;
		for (int i = 0; i < parent.noPix; i++) {
			pixnorm=0;
			for (int lambda = 0; lambda < parent.n; lambda++) {
				pixnorm+=parent.Xsub[i][lambda];
			}
			scalarprod2=0;scalarprod1=0;
			for (int lambda = 0; lambda < parent.n; lambda++) {
				scalarprod1+=(parent.Xsub[i][lambda]/pixnorm-base3[lambda])*base1[lambda];
				scalarprod2+=(parent.Xsub[i][lambda]/pixnorm-base3[lambda])*base2[lambda];
			}			
			s1=Ginv[0][0]*scalarprod1+Ginv[0][1]*scalarprod2;
			s2=Ginv[1][0]*scalarprod1+Ginv[1][1]*scalarprod2;
			i1=index_x(s1, s2);
			i2=index_y(s1, s2);
			if (i2>=0 && i2<maxpix && i1>=0 && i1<maxpix) conc_hist[i1][i2]++;
		}				
	}

	private void updateImage()
	{
		double maximage=0;
		for (int i = 0; i < conc_hist.length; i++) {
			for (int j = 0; j < conc_hist.length; j++) {
				conc_hist[i][j]=(float) Math.log(conc_hist[i][j]+1);
				if (conc_hist[i][j]>maximage) maximage=conc_hist[i][j];
			}
		}
		for (int i = 0; i < conc_hist.length; i++) {
			for (int j = 0; j < conc_hist.length; j++) {
				((byte [])triangle_image.getProcessor().getPixels())[i*maxpix+j]=((byte)(255*conc_hist[i][j]/maximage));
			}
		}
		for (int i = 0; i < conc_hist.length; i++) {
			for (int j = 0; j < conc_hist.length; j++) {
				if (spec_triangle[i][j]>0)
					((byte [])triangle_image.getProcessor().getPixels())[i*maxpix+j]=(byte)200;
			}
		}
		for (int i = 0; i < conc_hist.length; i++) {
			for (int j = 0; j < conc_hist.length; j++) {
				if (positivity_boundary[i][j]>0)
					((byte [])triangle_image.getProcessor().getPixels())[i*maxpix+j]+=(byte)50;
			}
		}
		if (canvas==null){
			triangle_image.createLut();
			triangle_image.show();
			canvas=triangle_image.getCanvas();
			canvas.addMouseListener(this);			
		}else{
			triangle_image.updateAndDraw();
		}
	}

	private int index_x(double s1, double s2) {
		return (int) (maxpix*((G[0][0]*s1+G[1][1]*Math.cos(alpha)*s2)/xscale-xoffset));
	}
	private int index_y(double s1, double s2) {
		return (int) (maxpix*((G[1][1]*Math.sin(alpha)*s2)/yscale-yoffset));
	}

	private double invindex_x(int i1, int i2) {
		double ii1=i1;
		double ii2=i2;
		return  ((ii1+xoffset*maxpix)*xscale*G[1][1]*Math.sin(alpha)-(ii2+yoffset*maxpix)*yscale*G[1][1]*Math.cos(alpha))/G[0][0]/G[1][1]/Math.sin(alpha)/maxpix;
	}
	private double invindex_y(int i1, int i2) {
		double ii2=i2;
		return (ii2+yoffset*maxpix)/G[1][1]/Math.sin(alpha)/maxpix*yscale;
	}

	private void calc_base() {
		for (int lambda = 0; lambda < parent.n; lambda++) {
			base1[lambda]=parent.modifiedS[dyes[0]][lambda]-parent.modifiedS[dyes[2]][lambda];
			base2[lambda]=parent.modifiedS[dyes[1]][lambda]-parent.modifiedS[dyes[2]][lambda];
			base3[lambda]=parent.modifiedS[dyes[2]][lambda];
		}
		G[0][0]=0;G[1][0]=0;G[0][1]=0;G[1][1]=0;
		for (int lambda = 0; lambda < parent.n; lambda++) {
			G[0][0]+=base1[lambda]*base1[lambda];
			G[0][1]+=base1[lambda]*base2[lambda];
			G[1][0]+=base2[lambda]*base1[lambda];
			G[1][1]+=base2[lambda]*base2[lambda];
		}
		double det=G[0][0]*G[1][1]-G[0][1]*G[1][0];
		Ginv[0][0]=G[1][1]/det;
		Ginv[1][1]=G[0][0]/det;
		Ginv[0][1]=-G[1][0]/det;
		Ginv[1][0]=-G[0][1]/det;
		alpha=Math.acos(G[0][1]/Math.sqrt(G[0][0]*G[1][1]));
	}

	private void default_values() {
		base1=new double [parent.n];
		base2=new double [parent.n];
		base3=new double [parent.n];
		G=new double [2][2];
		Ginv=new double [2][2];
	}

	//Mouse listener stuff
	public void mouseClicked(MouseEvent e) { 
	} 
	public void mousePressed(MouseEvent e) {
		int x = e.getX();
		int y = e.getY(); 
		int offscreenX = canvas.offScreenX(x);
		int offscreenY = canvas.offScreenY(y); 
		double s1=invindex_x(offscreenY,offscreenX);
		double s2=invindex_y(offscreenY,offscreenX);
		spectraChosen=chooseSpec(s1, s2);
		if (spectraChosen>=0) pressed=true;
		else pressed=false;
//		IJ.write("mousePressed: "+offscreenX+","+offscreenY+" "+s1+" "+s2); 
		
	} 
	public void mouseReleased(MouseEvent e) {
		if (pressed){
			int x = e.getX();
			int y = e.getY(); 
			int offscreenX = canvas.offScreenX(x);
			int offscreenY = canvas.offScreenY(y); 
			double s1=invindex_x(offscreenY,offscreenX);
			double s2=invindex_y(offscreenY,offscreenX);
			calc_new_spec(s1,s2);
			calc_and_draw();
			parent.plotSpectra();
			parent.showModifiedConcentrations();
//			IJ.write("mouse moved to : "+offscreenX+","+offscreenY+" "+x+" "+y+"  "+s1+" "+s2); 			
		}
	} 
	
	public void mouseEntered(MouseEvent e) {} 
	public void mouseExited(MouseEvent e) {} 

	private int chooseSpec(double s1, double s2){
		if (s1*s1+s2*s2<0.02) return 2;
		else if ((s1-1)*(s1-1)+s2*s2<0.02) return 0;
		else if ((s2-1)*(s2-1)+s1*s1<0.02) return 1;
		else return -1;
	}

	private void calc_new_spec(double s1, double s2){
		double [] s = new double [dyes.length];
		s[2]=1.0-s1-s2; s[0]=s1; s[1]=s2;
		double temp;
		for (int ch=0; ch<parent.n; ch++){
			temp=0;
			for (int dye=0; dye<dyes.length; dye++){
				temp+=parent.modifiedS[dyes[dye]][ch]*s[dye];
			}
			parent.modifiedS[dyes[spectraChosen]][ch]=temp;
		}
	}
}
