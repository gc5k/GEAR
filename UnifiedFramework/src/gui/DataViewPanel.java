package gui;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.Arrays;
import java.util.List;
import java.util.Vector;

import javax.swing.JScrollPane;
import javax.swing.JTable;
import javax.swing.table.DefaultTableModel;

public class DataViewPanel extends JScrollPane {

	//

	private DefaultTableModel tableModel;
	private JTable table;

	//

	DataViewPanel() {
		ini();
	}

	void ini() {
		tableModel = new DataViewTableModel();
		table = new JTable(tableModel);
		table.getTableHeader().setReorderingAllowed(false);
		setViewportView(table);
	}

	void setFile(File file) throws IOException {
		Vector data = new Vector();
		Vector cols = new Vector();
		BufferedReader br = new BufferedReader(new FileReader(file));
		String line;
		int index = 0;
		int colCount = 0;
		while ((line = br.readLine()) != null) {
			String[] cells = line.split("\\s");
			colCount = Math.max(colCount, cells.length);
			List<String> list = Arrays.asList(cells);
			Vector dataRow = new Vector(list);
			data.addElement(dataRow);
			if (++index == 100) {
				break;
			}
		}
		br.close();
		for (int i = 0; i < colCount; i++) {
			cols.addElement(String.valueOf(i + 1));
		}
		tableModel.setDataVector(data, cols);
	}

}

class DataViewTableModel extends DefaultTableModel {

	public boolean isCellEditable(int row, int column) {
		return false;
	}

}
