package gui;

import java.io.IOException;
import java.util.Vector;

import javax.swing.JScrollPane;
import javax.swing.JTable;
import javax.swing.table.DefaultTableModel;

public class DataViewPanel extends JScrollPane {

	//

	private JTable table;

	//

	DataViewPanel() {
		ini();
	}

	void ini() {
		table = new JTable();
		table.getTableHeader().setReorderingAllowed(false);
		setViewportView(table);
	}

	void setDataFile(DataFile dataFile) throws IOException {
		Vector[] vs = dataFile.getData(100);
		DefaultTableModel tableModel = new DataViewTableModel();
		tableModel.setDataVector(vs[0], vs[1]);
		table.setModel(tableModel);
	}

}

class DataViewTableModel extends DefaultTableModel {

	public boolean isCellEditable(int row, int column) {
		return false;
	}

}
