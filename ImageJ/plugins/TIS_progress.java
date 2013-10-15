import java.awt.*;
import java.awt.event.*;
import ij.plugin.frame.*;
import ij.*;
import ij.gui.*;

/**
*	Adapted from Image Processing Demo.
*/
public class TIS_progress extends PlugInFrame implements ActionListener {

	/**
	 * 
	 */
	private static final long serialVersionUID = 3728787123075205236L;
	Panel panel;
	int previousID;
	static Frame instance;
	TestImagingParameters_ TISinstance;
	
	public TIS_progress(TestImagingParameters_ TIS) {
		super("Modify channels");
		TISinstance=TIS;
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
		panel.setLayout(new GridLayout(1,2,2,2));
		for (int exc=0; exc<TIS.nexc; exc++)
		{
			addButton("Exc_"+Integer.toString(exc+1));
		}
		addButton("Ok");

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
		if (label=="Ok")
		{
			TISinstance.cancelled=true;
			dispose();
			instance = null;			
		}
		else
		{
			try
			{
				TISinstance.channelDialog(Integer.valueOf(label.substring(4))-1);
				TISinstance.calc_and_display();
			}
			finally
			{
			}
		}
	}

	public void processWindowEvent(WindowEvent e) {
		super.processWindowEvent(e);
		if (e.getID()==WindowEvent.WINDOW_CLOSING) {
			TISinstance.cancelled=true;
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
