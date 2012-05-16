package gui;

import java.awt.BorderLayout;
import java.awt.Container;
import java.awt.Dimension;
import java.awt.FlowLayout;
import java.awt.Font;
import java.awt.Frame;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.util.Enumeration;

import javax.swing.BorderFactory;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JDialog;
import javax.swing.JPanel;
import javax.swing.JTabbedPane;
import javax.swing.SwingUtilities;
import javax.swing.UIManager;
import javax.swing.plaf.FontUIResource;

public class DataLoadDialog extends JDialog {

	//

	JCheckBox checkBox;
	DataPanel dataPanel_ap;

	//

	public DataLoadDialog(Frame owner) {
		super(owner, true);
		ini();
	}

	private void ini() {
		//
		Container container = getContentPane();
		container.setLayout(new BorderLayout());
		//
		JTabbedPane tabbedPane = new JTabbedPane();
		tabbedPane.setBorder(BorderFactory.createRaisedBevelBorder());
		//
		DataPanel dataPanel_bi = new DataPanel();
		dataPanel_bi.setFileFormatList(new String[] { "bed", "bim", "fam" });
		dataPanel_bi.ini();
		tabbedPane.add("Binary Input", dataPanel_bi);
		//
		DataPanel dataPanel_si = new DataPanel();
		dataPanel_si.setFileFormatList(new String[] { "ped", "map" });
		dataPanel_si.ini();
		tabbedPane.add("Standard Input", dataPanel_si);
		//
		JPanel panel_ap = new JPanel();
		panel_ap.setLayout(new BorderLayout());
		JPanel panel_tmp = new JPanel();
		panel_tmp.setLayout(new GridBagLayout());
		checkBox = new JCheckBox();
		checkBox.setText("Use Alternative Phenotype");
		checkBox.addItemListener(new ItemListener() {
			public void itemStateChanged(ItemEvent e) {
				DataLoadDialog.this.checkBox_itemStateChanged(e);
			}
		});
		panel_tmp.add(checkBox, new GridBagConstraints(0, 0, 1, 1, 0.0, 1.0, GridBagConstraints.CENTER, GridBagConstraints.BOTH, new Insets(5, 5, 5, 5), 0, 0));
		panel_tmp.add(new JPanel(), new GridBagConstraints(1, 0, 1, 1, 1.0, 1.0, GridBagConstraints.CENTER, GridBagConstraints.BOTH, new Insets(5, 5, 5, 5), 0,
				0));
		panel_ap.add(panel_tmp, BorderLayout.NORTH);
		//
		dataPanel_ap = new DataPanel();
		dataPanel_ap.setFileFormatList(new String[] { "phe" });
		dataPanel_ap.ini();
		panel_ap.add(dataPanel_ap, BorderLayout.CENTER);
		//
		tabbedPane.add("Alternate Phenotype", panel_ap);
		//
		container.add(tabbedPane, BorderLayout.NORTH);
		//
		JPanel panel_c = new JPanel();
		panel_c.setLayout(new FlowLayout(FlowLayout.CENTER, 5, 25));
		JButton button = new JButton();
		button.setText("OK");
		button.setPreferredSize(new Dimension(120, 30));
		button.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				DataLoadDialog.this.dispose();
			}
		});
		panel_c.add(button);
		container.add(panel_c, BorderLayout.CENTER);
		//
		setTitle("Load Data");
		setResizable(false);
		setDefaultCloseOperation(DISPOSE_ON_CLOSE);
		//
		checkBox_itemStateChanged(null);
	}

	public void display() {
		setSize(800, 400);
		setLocationRelativeTo(getOwner());
		setVisible(true);
	}

	private void checkBox_itemStateChanged(ItemEvent e) {
		dataPanel_ap.setEnabled_2(checkBox.isSelected());
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
				new DataLoadDialog(null).display();
			}
		});
	}

}
