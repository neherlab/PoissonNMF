import java.awt.*;
import java.awt.event.*;
import ij.plugin.frame.*;
import ij.*;
import ij.gui.*;

/**
*	Adapted from Image Processing Demo.
*/
public class OptimizeAcquisitionProgress extends PlugInFrame implements ActionListener {

	/**
	 * 
	 */
	private static final long serialVersionUID = 3728787123075205236L;
	Panel panel;
	int previousID;
	static Frame instance;
	OptimizeAcquisitionScheme_ OASinstance;
	
	public OptimizeAcquisitionProgress(OptimizeAcquisitionScheme_ OAS) {
		super("OAS is running...");
		OASinstance=OAS;
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
		panel.setLayout(new GridLayout(1, 1, 2, 2));
		addButton("Cancel");
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
		if (label=="Cancel")
		{
			OASinstance.cancelled=true;
			dispose();
			instance = null;			
		}
	}

	public void processWindowEvent(WindowEvent e) {
		super.processWindowEvent(e);
		if (e.getID()==WindowEvent.WINDOW_CLOSING) {
			OASinstance.cancelled=true;
			dispose();
			instance = null;			
		}
	}	
	public void closeWindow()
	{
		dispose();
		instance = null;					
	}
	
	
}
