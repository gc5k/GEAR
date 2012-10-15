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
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JTabbedPane;
import javax.swing.SwingUtilities;
import javax.swing.UIManager;
import javax.swing.plaf.FontUIResource;

import org.apache.commons.lang3.ArrayUtils;

public class DataLoadDialog extends JDialog {

	//

	private DataLoadPanel dataPanel_bi;
	private DataLoadPanel dataPanel_si;
	private DataLoadPanel dataPanel_ap;
	private JCheckBox useAP_checkBox;

	private DataFile[][] dataFiles = new DataFile[2][];

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
		dataPanel_bi = new DataLoadPanel();
		dataPanel_bi.setFileFormats(new String[] { "bed", "bim", "fam" });
		dataPanel_bi.ini();
		tabbedPane.add("Binary Input", dataPanel_bi);
		//
		dataPanel_si = new DataLoadPanel();
		dataPanel_si.setFileFormats(new String[] { "ped", "map" });
		dataPanel_si.ini();
		tabbedPane.add("Standard Input", dataPanel_si);
		//
		JPanel panel_ap = new JPanel();
		panel_ap.setLayout(new BorderLayout());
		JPanel panel_tmp = new JPanel();
		panel_tmp.setLayout(new GridBagLayout());
		useAP_checkBox = new JCheckBox();
		useAP_checkBox.setText("Use Alternative Phenotype");
		useAP_checkBox.addItemListener(new ItemListener() {
			public void itemStateChanged(ItemEvent e) {
				DataLoadDialog.this.checkBox_itemStateChanged(e);
			}
		});
		panel_tmp.add(useAP_checkBox, new GridBagConstraints(0, 0, 1, 1, 0.0, 1.0, GridBagConstraints.CENTER, GridBagConstraints.BOTH, new Insets(10, 5, 5, 5),
				0, 0));
		panel_ap.add(panel_tmp, BorderLayout.NORTH);
		//
		dataPanel_ap = new DataLoadPanel();
		dataPanel_ap.setFileFormats(new String[] { "phe" });
		dataPanel_ap.ini();
		panel_ap.add(dataPanel_ap, BorderLayout.CENTER);
		//
		tabbedPane.add("Alternate Phenotype", panel_ap);
		//
		container.add(tabbedPane, BorderLayout.NORTH);
		//
		JPanel panel_c = new JPanel();
		panel_c.setLayout(new FlowLayout(FlowLayout.CENTER, 5, 30));
		JButton button = new JButton();
		button.setText("OK");
		button.setPreferredSize(new Dimension(120, 30));
		button.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				DataLoadDialog.this.okFiles();
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
		dataPanel_ap.setEnabled_2(useAP_checkBox.isSelected());
	}

	private void okFiles() {
		DataFile[] files_bi = dataPanel_bi.getDataFiles();
		if (!ArrayUtils.isEmpty(files_bi)) {
			dataFiles[0] = files_bi;
		}
		DataFile[] files_si = dataPanel_si.getDataFiles();
		if (!ArrayUtils.isEmpty(files_si)) {
			dataFiles[0] = files_si;
		}
		DataFile[] files_ap = dataPanel_ap.getDataFiles();
		if (!ArrayUtils.isEmpty(files_ap)) {
			dataFiles[1] = files_ap;
		}
		if (ArrayUtils.isEmpty(dataFiles[0]) && ArrayUtils.isEmpty(dataFiles[1])) {
			JOptionPane.showMessageDialog(this, "No data input!", "Error", JOptionPane.ERROR_MESSAGE);
			return;
		}
		dispose();
	}

	public DataFile[][] getDataFiles() {
		return dataFiles;
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
