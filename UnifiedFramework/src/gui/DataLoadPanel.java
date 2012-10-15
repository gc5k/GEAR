package gui;

import java.awt.Dimension;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.io.File;
import java.io.FilenameFilter;
import java.util.TreeSet;

import javax.swing.BorderFactory;
import javax.swing.DefaultComboBoxModel;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JComboBox;
import javax.swing.JFileChooser;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JTextField;

public class DataLoadPanel extends JPanel {

	//

	private DataFile[] dataFiles;
	private JTextField[] fileTextFieldList;
	private JCheckBox[] fileCheckBoxList;
	private JButton[] fileButtonList;

	private File dir;
	private JCheckBox dirCheckBox;
	private JTextField dirTextField;
	private JButton dirButton;
	private JComboBox dirComboBox;

	//

	public void setFileFormats(String[] formatList) {
		dataFiles = new DataFile[formatList.length];
		for (int i = 0; i < dataFiles.length; i++) {
			dataFiles[i] = new DataFile();
			dataFiles[i].setFormat(formatList[i]);
		}
		fileTextFieldList = new JTextField[formatList.length];
		fileCheckBoxList = new JCheckBox[formatList.length];
		fileButtonList = new JButton[formatList.length];
	}

	public void ini() {
		setLayout(new GridBagLayout());
		//
		JPanel panel_1 = new JPanel();
		panel_1.setLayout(new GridBagLayout());
		for (int i = 0; i < dataFiles.length; i++) {
			panel_1.add(new JLabel("." + dataFiles[i].getFormat() + " file"), new GridBagConstraints(0, i, 1, 1, 0.0, 1.0, GridBagConstraints.CENTER,
					GridBagConstraints.BOTH, new Insets(5, 5, 5, 5), 0, 0));
			fileTextFieldList[i] = new JTextField();
			fileTextFieldList[i].setEditable(false);
			panel_1.add(fileTextFieldList[i], new GridBagConstraints(1, i, 1, 1, 1.0, 1.0, GridBagConstraints.CENTER, GridBagConstraints.HORIZONTAL,
					new Insets(5, 5, 5, 5), 0, 0));
			fileCheckBoxList[i] = new JCheckBox();
			fileCheckBoxList[i].setText("has header");
			panel_1.add(fileCheckBoxList[i], new GridBagConstraints(2, i, 1, 1, 0.0, 1.0, GridBagConstraints.CENTER, GridBagConstraints.HORIZONTAL, new Insets(
					5, 5, 5, 5), 0, 0));
			fileButtonList[i] = new JButton();
			fileButtonList[i].setText("Browse");
			final int j = i;
			fileButtonList[i].addActionListener(new ActionListener() {
				public void actionPerformed(ActionEvent e) {
					DataLoadPanel.this.button_actionPerformed(e, j);
				}
			});
			panel_1.add(fileButtonList[i], new GridBagConstraints(3, i, 1, 1, 0.0, 1.0, GridBagConstraints.CENTER, GridBagConstraints.HORIZONTAL, new Insets(5,
					5, 5, 5), 0, 0));
		}
		add(panel_1, new GridBagConstraints(0, 0, 1, 1, 1.0, 0.0, GridBagConstraints.CENTER, GridBagConstraints.BOTH, new Insets(5, 5, 5, 5), 0, 0));
		//
		JPanel panel_qf = new JPanel(); // quick fileset
		panel_qf.setBorder(BorderFactory.createTitledBorder("Quick Fileset"));
		panel_qf.setLayout(new GridBagLayout());
		dirCheckBox = new JCheckBox();
		dirCheckBox.addItemListener(new ItemListener() {
			public void itemStateChanged(ItemEvent e) {
				DataLoadPanel.this.dirCheckBox_itemStateChanged(e);
			}
		});
		dirCheckBox.setSelected(false);
		panel_qf.add(dirCheckBox,
				new GridBagConstraints(0, 0, 1, 1, 0.0, 1.0, GridBagConstraints.CENTER, GridBagConstraints.BOTH, new Insets(5, 5, 5, 5), 0, 0));
		dirTextField = new JTextField();
		dirTextField.setEditable(false);
		panel_qf.add(dirTextField, new GridBagConstraints(1, 0, 1, 1, 1.0, 1.0, GridBagConstraints.CENTER, GridBagConstraints.BOTH, new Insets(5, 5, 5, 5), 0,
				0));
		dirButton = new JButton();
		dirButton.setText("...");
		dirButton.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				DataLoadPanel.this.button_actionPerformed(e, dataFiles.length);
			}
		});
		panel_qf.add(dirButton, new GridBagConstraints(2, 0, 1, 1, 0.0, 1.0, GridBagConstraints.CENTER, GridBagConstraints.BOTH, new Insets(5, 5, 5, 5), 0, 0));
		dirComboBox = new JComboBox();
		dirComboBox.setPreferredSize(new Dimension(150, 28));
		dirComboBox.setEditable(false);
		dirComboBox.addItemListener(new ItemListener() {
			public void itemStateChanged(ItemEvent e) {
				DataLoadPanel.this.dirComboBox_itemStateChanged(e);
			}
		});
		panel_qf.add(dirComboBox,
				new GridBagConstraints(3, 0, 1, 1, 0.0, 1.0, GridBagConstraints.CENTER, GridBagConstraints.BOTH, new Insets(5, 5, 5, 5), 0, 0));
		add(panel_qf, new GridBagConstraints(0, 1, 1, 1, 1.0, 0.0, GridBagConstraints.CENTER, GridBagConstraints.BOTH, new Insets(5, 5, 5, 5), 0, 0));
		//
		add(new JPanel(), new GridBagConstraints(0, 2, 1, 1, 1.0, 1.0, GridBagConstraints.CENTER, GridBagConstraints.BOTH, new Insets(5, 5, 5, 5), 0, 0));
		//
		dirCheckBox_itemStateChanged(null);
	}

	private void dirComboBox_itemStateChanged(ItemEvent e) {
		if (e != null && e.getStateChange() == ItemEvent.DESELECTED) {
			return;
		}
		clearFiles();
		final String fileprefix = (String) dirComboBox.getSelectedItem();
		File[] files = dir.listFiles(new FilenameFilter() {
			public boolean accept(File dir, String fileName) {
				return fileName.startsWith(fileprefix + ".");
			}
		});
		for (int i = 0; i < dataFiles.length; i++) {
			for (int j = 0; j < files.length; j++) {
				if (files[j].getName().toLowerCase().endsWith("." + dataFiles[i].getFormat().toLowerCase())) {
					setFile(files[j], i);
				}
			}
		}
	}

	private void button_actionPerformed(ActionEvent e, int id) {
		if (-1 < id && id < dataFiles.length) {
			JFileChooser fileChooser = new JFileChooser();
			fileChooser.setFileSelectionMode(JFileChooser.FILES_ONLY);
			fileChooser.setMultiSelectionEnabled(false);
			if (fileChooser.showOpenDialog(this) == JFileChooser.APPROVE_OPTION) {
				setFile(fileChooser.getSelectedFile(), id);
			}
			return;
		}
		if (id == dataFiles.length) {
			JFileChooser fileChooser = new JFileChooser();
			fileChooser.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
			fileChooser.setMultiSelectionEnabled(false);
			if (fileChooser.showOpenDialog(this) == JFileChooser.APPROVE_OPTION) {
				setDir(fileChooser.getSelectedFile());
				File[] files = dir.listFiles(new FilenameFilter() {
					public boolean accept(File dir, String fileName) {
						for (DataFile dataFile : dataFiles) {
							if (fileName.toLowerCase().endsWith("." + dataFile.getFormat().toLowerCase())) {
								return true;
							}
						}
						return false;
					}
				});
				TreeSet<String> fileprefix1 = new TreeSet<String>();
				for (File file : files) {
					String fileName = file.getName();
					fileName = fileName.substring(0, fileName.lastIndexOf('.'));
					fileprefix1.add(fileName);
				}
				String[] fileprefix2 = fileprefix1.toArray(new String[fileprefix1.size()]);
				dirComboBox.setModel(new DefaultComboBoxModel(fileprefix2));
				dirComboBox_itemStateChanged(null);
			}
			return;
		}
	}

	private void dirCheckBox_itemStateChanged(ItemEvent e) {
		for (int i = 0; i < fileButtonList.length; i++) {
			fileButtonList[i].setEnabled(!dirCheckBox.isSelected());
		}
		setEnabled_1(dirCheckBox.isSelected());
	}

	private void setFile(File file, int id) {
		if (file == null) {
			dataFiles[id].setFile(null);
			fileTextFieldList[id].setText(null);
			fileTextFieldList[id].setToolTipText(null);
		} else {
			dataFiles[id].setFile(file);
			fileTextFieldList[id].setText(file.getAbsolutePath());
			fileTextFieldList[id].setToolTipText(file.getAbsolutePath());
		}
	}

	private void clearFiles() {
		for (int i = 0; i < dataFiles.length; i++) {
			setFile(null, i);
		}
	}

	private void setDir(File dir) {
		this.dir = dir;
		dirTextField.setText(dir.getAbsolutePath());
		dirTextField.setToolTipText(dir.getAbsolutePath());
	}

	private void setEnabled_1(boolean enabled) {
		dirTextField.setEnabled(enabled);
		dirButton.setEnabled(enabled);
		dirComboBox.setEnabled(enabled);
	}

	void setEnabled_2(boolean enabled) {
		for (int i = 0; i < dataFiles.length; i++) {
			fileTextFieldList[i].setEnabled(enabled);
			fileCheckBoxList[i].setEnabled(enabled);
		}
		dirCheckBox.setEnabled(enabled);
		if (enabled) {
			dirCheckBox_itemStateChanged(null);
		} else {
			for (int i = 0; i < dataFiles.length; i++) {
				fileButtonList[i].setEnabled(false);
			}
			setEnabled_1(false);
		}
	}

	public DataFile[] getDataFiles() {
		return dataFiles;
	}

}
