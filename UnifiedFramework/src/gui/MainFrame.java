package gui;

import java.awt.Font;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.Enumeration;

import javax.swing.JFrame;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;
import javax.swing.JPanel;
import javax.swing.JTabbedPane;
import javax.swing.SwingUtilities;
import javax.swing.UIManager;
import javax.swing.plaf.FontUIResource;

public class MainFrame extends JFrame {

	//

	private DataLoadDialog dataLoadFrame;

	//

	MainFrame() {
		super();
		// JMenuBar
		JMenuBar menubar = new JMenuBar();

		JMenu menu_project = new JMenu("Project");
		menubar.add(menu_project);
		JMenuItem menuItem_np = new JMenuItem("New");
		menu_project.add(menuItem_np);
		JMenuItem menuItem_op = new JMenuItem("Open");
		menu_project.add(menuItem_op);
		JMenuItem menuItem_sp = new JMenuItem("Save");
		menu_project.add(menuItem_sp);
		JMenuItem menuItem_rp = new JMenuItem("Save As");
		menu_project.add(menuItem_rp);
		menu_project.addSeparator();
		JMenuItem menuItem_acon = new JMenuItem("Configurations");
		menu_project.add(menuItem_acon);
		menu_project.addSeparator();
		JMenuItem menuItem_recent = new JMenuItem("Recent Project(s)");
		menu_project.add(menuItem_recent);
		menu_project.addSeparator();
		JMenuItem menuItem_exit = new JMenuItem("Exit");
		menuItem_exit.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				MainFrame.this.actionPerformed_exit();
			}
		});
		menu_project.add(menuItem_exit);

		JMenu menu_data = new JMenu("Data");
		menubar.add(menu_data);
		JMenuItem menuItem_ld = new JMenuItem("Load");
		menuItem_ld.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				MainFrame.this.actionPerformed_dataload();
			}
		});
		menu_data.add(menuItem_ld);
		JMenuItem menuItem_od = new JMenuItem("Output");
		menuItem_od.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				MainFrame.this.actionPerformed_dataoutput();
			}
		});
		menu_data.add(menuItem_od);
		JMenuItem menuItem_md = new JMenuItem("Merge");
		menu_data.add(menuItem_md);
		menu_data.addSeparator();
		JMenuItem menuItem_vd = new JMenuItem("View");
		menu_data.add(menuItem_vd);
		JMenuItem menuItem_dds = new JMenuItem("Define Data Format");
		menu_data.add(menuItem_dds);
		JMenuItem menuItem_fd = new JMenuItem("Filter");
		menu_data.add(menuItem_fd);

		JMenu menu_tools = new JMenu("Tools");
		menubar.add(menu_tools);
		JMenuItem menuItem_ss = new JMenuItem("Summary Statistics");
		menu_tools.add(menuItem_ss);

		JMenu menu_Advanced = new JMenu("Advanced");
		menubar.add(menu_Advanced);
		JMenuItem menuItem_igo = new JMenuItem("Import GMDR Operation");
		menu_Advanced.add(menuItem_igo);
		JMenuItem menuItem_cgo = new JMenuItem("Create GMDR Operation");
		menu_Advanced.add(menuItem_cgo);
		JMenuItem menuItem_ano = new JMenuItem("Add non-GMDR Operation");
		menu_Advanced.add(menuItem_ano);

		JMenu menu_Help = new JMenu("Help");
		menubar.add(menu_Help);
		JMenuItem menuItem_a = new JMenuItem("About");
		menu_Help.add(menuItem_a);
		JMenuItem menuItem_h = new JMenuItem("Help");
		menu_Help.add(menuItem_h);

		setJMenuBar(menubar);
		//
		JTabbedPane tabbedPane = new JTabbedPane();
		tabbedPane.add("MDR Analysis", new JPanel());
		tabbedPane.add("Metric Calculation", new JPanel());
		tabbedPane.add("Settings", new JPanel());
		setContentPane(tabbedPane);
		//
		setDefaultCloseOperation(DISPOSE_ON_CLOSE);
		setTitle("Generalized Multifactor Dimensionality Reduction V1.0");
		setSize(800, 600);
		setLocationRelativeTo(null);
		setExtendedState(MAXIMIZED_BOTH);
	}

	public void dispose() {
		super.dispose();
		System.exit(0);
	}

	private void actionPerformed_dataload() {
		if (dataLoadFrame == null) {
			dataLoadFrame = new DataLoadDialog(this);
		}
		dataLoadFrame.display();
	}

	private void actionPerformed_dataoutput() {
		if (dataLoadFrame == null) {
			dataLoadFrame = new DataLoadDialog(this);
		}
		dataLoadFrame.display();
	}

	private void actionPerformed_exit() {
		dispose();
	}

	public static void main(String[] args) {
		SwingUtilities.invokeLater(new Runnable() {
			public void run() {
				try {
					// UIManager.setLookAndFeel(UIManager.getCrossPlatformLookAndFeelClassName());
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
				new MainFrame().setVisible(true);
			}
		});
	}

}
