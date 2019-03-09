package gear.family.plink;

import gear.family.pedigree.file.MapFile;
import gear.family.pedigree.file.PedigreeFile;

public class PLINKTextParser extends PLINKParser {
	public PLINKTextParser(String ped, String map) {
		super(new PedigreeFile(ped), new MapFile(map));
	}

	@Override
	public void parse() {
		super.parse();
		mapData.setPolymorphismMarker(pedigreeData.getPolymorphism());
		if (mapData.getFilename() == null) {
			mapData.setMarker(pedigreeData.getNumMarkers());
		}
	}
}
