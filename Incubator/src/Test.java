import java.io.FileWriter;
import java.io.BufferedWriter;
import java.text.ParseException;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.PosixParser;

public class Test {

	public static void main(String[] args) throws org.apache.commons.cli.ParseException {

		CommandLineParser parser = new PosixParser();

		// create the Options
		Options options = new Options();
		options.addOption( "a", "all", true, "do not hide entries starting with ." );
		options.addOption( "A", "almost-all", false, "do not list implied . and .." );
		options.addOption( "b", "escape", false, "print octal escapes for nongraphic "
		                                         + "characters" );
		options.addOption( OptionBuilder.withLongOpt( "block-size" )
		                                .withDescription( "use SIZE-byte blocks" )
		                                .hasArg()
		                                .withArgName("SIZE")
		                                .create('k') );
		options.addOption( "B", "ignore-backups", false, "do not list implied entried "
		                                                 + "ending with ~");
		options.addOption( "c", false, "with -lt: sort by, and show, ctime (time of last " 
		                               + "modification of file status information) with "
		                               + "-l:show ctime and sort by name otherwise: sort "
		                               + "by ctime" );
		options.addOption( "C", false, "list entries by columns" );

		String[] arg = new String[]{ "--all=3 ", "--block-size=10 5"  };

		CommandLine line = null;
		line = parser.parse( options, args );
		if( line.hasOption( "block-size" ) ) {
				// print the value of block-size
			String[] v = line.getOptionValues("block-size");
			for(int i = 0; i < v.length; i++) {
				System.out.println(v[i]);
			}
//				System.out.println( line.getOptionValues( "block-size" ) );
		}
		if( line.hasOption("a")) {
			System.out.println(line.getOptionValue("a"));
		}
	}
}
