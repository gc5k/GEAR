package gui;

import java.awt.BorderLayout;
import java.awt.Container;
import java.awt.Font;
import java.awt.Frame;
import java.io.File;
import java.io.IOException;
import java.util.Enumeration;

import javax.swing.BorderFactory;
import javax.swing.JDialog;
import javax.swing.JOptionPane;
import javax.swing.JTabbedPane;
import javax.swing.SwingUtilities;
import javax.swing.UIManager;
import javax.swing.plaf.FontUIResource;

public class DataViewDialog extends JDialog {

	//

	private File[] files;

	//

	public DataViewDialog(Frame owner) {
		super(owner, true);
	}

	void ini(File[] files) {
		this.files = files;
		//
		Container container = getContentPane();
		container.setLayout(new BorderLayout());
		//
		JTabbedPane tabbedPane = new JTabbedPane();
		tabbedPane.setBorder(BorderFactory.createRaisedBevelBorder());
		//
		for (File file : this.files) {
			DataViewPanel dataViewPanel = new DataViewPanel();
			try {
				dataViewPanel.setFile(file);
			} catch (IOException e) {
				JOptionPane.showMessageDialog(this, e);
			}
			tabbedPane.add(file.getName(), dataViewPanel);
		}
		//
		container.add(tabbedPane, BorderLayout.CENTER);
		//
		setTitle("View Data (Top 100 rows)");
		setResizable(true);
		setDefaultCloseOperation(DISPOSE_ON_CLOSE);
	}

	public void display() {
		setSize(800, 600);
		setLocationRelativeTo(getOwner());
		setVisible(true);
	}

	public static void main(String[] args) {
		SwingUtilities.invokeLater(new Runnable() {
			public void run() {
				try {
					UIManager.setLookAndFeel("com.sun.java.swing.plaf.nimbus.NimbusLookAndFeel");
				} catch (Exception e) {
				}
				FontUIResource fontRes = new FontUIResource(new Font(Font.DIALOG, Font.PLAIN, 12));
				for (Enumeration<Object> keys = UIManager.getDefaults().keys(); keys.hasMoreElements();) {
					Object key = keys.nextElement();
					Object value = UIManager.get(key);
					if (value instanceof FontUIResource) {
						UIManager.put(key, fontRes);
					}
				}
				new DataViewDialog(null).display();
			}
		});
	}

}
