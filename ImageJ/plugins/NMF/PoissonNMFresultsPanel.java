import java.awt.*;
import java.awt.event.*;
import ij.plugin.frame.*;
import ij.*;
import ij.gui.*;

/**
*	Adapted from Image Processing Demo.
*/
public class PoissonNMFresultsPanel extends PlugInFrame implements ActionListener {

	/**
	 * 
	 */
	private static final long serialVersionUID = -5239603611196213480L;
	Panel panel;
	int previousID;
	static Frame instance;
	PoissonNMF_ pnmfinstance;
	
	public PoissonNMFresultsPanel(PoissonNMF_ pnmf) {
		super("PoissonNMF Results");
		pnmfinstance=pnmf;
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
		addButton("Save Spectra");
		addButton("Background Map");
		addButton("Background Spectrum");
		addButton("Simplex projection");
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
			pnmfinstance.RGB_overlay();
		else if  (label=="Save Spectra")
			pnmfinstance.saveSpectra();
		else if  (label=="Background Map")
			pnmfinstance.BG_map();
		else if  (label=="Background Spectrum")
			pnmfinstance.plotbg();		
		else if  (label=="Simplex projection")
			pnmfinstance.show_simplex();
		else if (label=="Exit")
		{
			dispose();
			pnmfinstance.plotw.close();
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
