import java.awt.*;
import java.awt.event.*;
import ij.plugin.frame.*;
import ij.*;
import ij.gui.*;

/**
*	Adapted from Image Processing Demo.
*/
public class PoissonNTFresultsPanel extends PlugInFrame implements ActionListener {

	/**
	 * 
	 */
	private static final long serialVersionUID = -5239603611196213480L;
	Panel panel;
	int previousID;
	static Frame instance;
	PoissonNTF_ pntfinstance;
	
	public PoissonNTFresultsPanel(PoissonNTF_ pntf) {
		super("PoissonNTF Results");
		pntfinstance=pntf;
		if (IJ.versionLessThan("1.39t"))
			return;
		if (instance!=null) {
			instance.toFront();
			return;
		}
		instance = this;
		addKeyListener(IJ.getInstance());
		setLayout(new FlowLayout());
		panel = new Panel();
		panel.setLayout(new GridLayout(3, 3, 2, 2));
		addButton("RGB overlay");
		addButton("Save Emission Spectra");
		addButton("Save Excitation Spectra");
		addButton("Background Map");
		addButton("Background Spectrum");
		addButton("Toggle Normalization");
		addButton("Exit");
		add(panel);

		pack();
		GUI.center(this);
		setVisible(true);
	}
	
	void addButton(String label) {
		Button b = new Button(label);
		b.addActionListener(this);
		b.addKeyListener(IJ.getInstance());
		panel.add(b);
	}

	public void actionPerformed(ActionEvent e) {
		String label = e.getActionCommand();
		if (label==null)
			return;
		if  (label=="RGB overlay")
			pntfinstance.RGB_overlay();
		else if  (label=="Save Emission Spectra")
			pntfinstance.saveEmissionSpectra();
		else if  (label=="Save Excitation Spectra")
			pntfinstance.saveExcitationSpectra();
		else if  (label=="Background Map")
			pntfinstance.BG_map();
		else if  (label=="Background Spectrum")
			pntfinstance.plotbg();		
		else if  (label=="Toggle Normalization")
			pntfinstance.togglenormalization();		
		else if (label=="Exit")
		{
			dispose();
			pntfinstance.plotEmission.close();
			pntfinstance.plotExcitation.close();
			instance=null;
		}
	}

	public void processWindowEvent(WindowEvent e) {
		super.processWindowEvent(e);
		if (e.getID()==WindowEvent.WINDOW_CLOSING) {
			dispose();
			instance = null;
			
		}
	}

}
