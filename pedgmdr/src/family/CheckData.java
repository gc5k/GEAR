package family;

/*
 * $Id: CheckData.java,v 3.14 2006/05/24 19:37:37 djbender Exp $
 * WHITEHEAD INSTITUTE
 * SOFTWARE COPYRIGHT NOTICE AGREEMENT
 * This software and its documentation are copyright 2003 by the
 * Whitehead Institute for Biomedical Research.  All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support
 * whatsoever.  The Whitehead Institute can not be responsible for its
 * use, misuse, or functionality.
 */
import java.util.ArrayList;
import java.util.Enumeration;
import java.util.Hashtable;
import java.util.Vector;

import edu.mit.wi.haploview.Chromosome;
import edu.mit.wi.pedfile.Family;
import edu.mit.wi.pedfile.Individual;
import edu.mit.wi.pedfile.MarkerResult;
import edu.mit.wi.pedfile.MendelError;
import edu.mit.wi.pedfile.PedFileException;

/**
 * <p>
 * Title: CheckData.java
 * </p>
 * <p>
 * Description: Used to check pedigree file and return information about observed heterozyosity, predicted
 * heterozygosity, Hardy-Weinberg test p-value, genotyped percent, number of families with a fully genotyped trio and
 * number of Mendelian inheritance errors.
 * </p>
 */
public class CheckData {

    private MDRPed pedFile;
    static public double hwCut = 0.001;
    static public int failedGenoCut = 75;
    static public int numMendErrCut = 1;
    static public double mafCut = 0.001;
    // These are unchanging defaults used for the reset button.
    static public double defaultHwCut = 0.001;
    static public int defaultFailedGenoCut = 75;
    static public int defaultNumMendErrCut = 1;
    static public double defaultMafCut = 0.001;
    static public int numAllele = 4;

    public CheckData(MDRPed pedFile) {
        this.pedFile = pedFile;
    }

    /**
     * checks the pedigree file
     * 
     * @return ArrayList includes information about observed heterozyosity, predicted heterozygosity, Hardy-Weinberg test
     *         p-value, genotyped percent, number of families with a fully genotyped trio and number of Mendelian
     *         inheritance errors.
     */
    public ArrayList check() throws MDRPedFileException, PedFileException {
        int numOfMarkers = pedFile.getNumMarkers();
        ArrayList results = new ArrayList(numOfMarkers);
        // _size = pedFile.getNumIndividuals();

        for (int i = 0; i < numOfMarkers; i++) {
            results.add(checkMarker(i));
        }
        return results;
    }

    private MarkerResult checkMarker(int loc) throws MDRPedFileException,
            PedFileException {
        MarkerResult result = new MarkerResult();
        Individual currentInd;
        int missing = 0, founderHetCount = 0, mendErrNum = 0;
        int allele1 = 0, allele2 = 0, hom = 0, het = 0, haploid = 0;
        Hashtable founderGenoCount = new Hashtable();
        Hashtable kidgeno = new Hashtable();
        int[] founderHomCount = new int[5];
        Vector mendels = new Vector();

        int[] count = new int[5];
        for (int i = 0; i < 5; i++) {
            founderHomCount[i] = 0;
            count[i] = 0;
        }

        // loop through each family, check data for marker loc
        Enumeration famList = pedFile.getFamList();
        while (famList.hasMoreElements()) {
            Family currentFamily = pedFile.getFamily((String) famList.nextElement());
            Enumeration indList = currentFamily.getMemberList();
            // loop through each individual in the current Family
            while (indList.hasMoreElements()) {
                currentInd = currentFamily.getMember((String) indList.nextElement());
                allele1 = currentInd.getAllele(loc, 0);
                allele2 = currentInd.getAllele(loc, 1);

                // if haploid, check for male hets
                if (Chromosome.getDataChrom().equals("chrx") && currentInd.getGender() == 1) {
                    if (allele1 != allele2) {
                        currentInd.zeroOutMarker(loc);
                        pedFile.setHaploidHets(true);
                    }
                }

                // no allele data missing
                if (allele1 > 0 && allele2 > 0) {
                    // make sure entry has parents
                    if (currentFamily.containsMember(currentInd.getMomID()) && currentFamily.containsMember(currentInd.getDadID())) {
                        // do mendel check
                        int momAllele1 = (currentFamily.getMember(currentInd.getMomID())).getAllele(loc, 0);
                        int momAllele2 = (currentFamily.getMember(currentInd.getMomID())).getAllele(loc, 1);
                        int dadAllele1 = (currentFamily.getMember(currentInd.getDadID())).getAllele(loc, 0);
                        int dadAllele2 = (currentFamily.getMember(currentInd.getDadID())).getAllele(loc, 1);

                        if (Chromosome.getDataChrom().equals("chrx")) {
                            if (currentInd.getGender() == 1) {
                                if (!(momAllele1 == 0 || momAllele2 == 0 || dadAllele1 == 0)) {
                                    // this is an x chrom for a male, so the only thing we need to check is if
                                    // allele1 matches either momallele1 or momallele2
                                    if (allele1 != momAllele1 && allele1 != momAllele2) {
                                        mendErrNum++;
                                        MendelError mend = new MendelError(
                                                currentInd.getFamilyID(),
                                                currentInd.getIndividualID());
                                        mendels.add(mend);
                                        currentInd.zeroOutMarker(loc);
                                        currentFamily.getMember(
                                                currentInd.getMomID()).zeroOutMarker(loc);
                                        currentFamily.getMember(
                                                currentInd.getDadID()).zeroOutMarker(loc);
                                    }
                                }
                            } else {
                                // if gender is anything except 1 we assume female
                                if (momAllele1 == momAllele2) {
                                    // mom hom and dad matches mom
                                    if (dadAllele1 == momAllele1) {
                                        // kid must be hom same allele
                                        if (allele1 != momAllele1 || allele2 != momAllele2) {
                                            mendErrNum++;
                                            MendelError mend = new MendelError(
                                                    currentInd.getFamilyID(),
                                                    currentInd.getIndividualID());
                                            mendels.add(mend);
                                            currentInd.zeroOutMarker(loc);
                                            currentFamily.getMember(
                                                    currentInd.getMomID()).zeroOutMarker(loc);
                                            currentFamily.getMember(
                                                    currentInd.getDadID()).zeroOutMarker(loc);
                                        }
                                    } else {
                                        // kid must be het
                                        if (allele1 == allele2) {
                                            mendErrNum++;
                                            MendelError mend = new MendelError(
                                                    currentInd.getFamilyID(),
                                                    currentInd.getIndividualID());
                                            mendels.add(mend);
                                            currentInd.zeroOutMarker(loc);
                                            currentFamily.getMember(
                                                    currentInd.getMomID()).zeroOutMarker(loc);
                                            currentFamily.getMember(
                                                    currentInd.getDadID()).zeroOutMarker(loc);
                                        }
                                    }
                                } else {
                                    // mom het,so only need to check that at least one allele matches dad
                                    if (allele1 != dadAllele1 && allele2 != dadAllele2) {
                                        mendErrNum++;
                                        MendelError mend = new MendelError(
                                                currentInd.getFamilyID(),
                                                currentInd.getIndividualID());
                                        mendels.add(mend);
                                        currentInd.zeroOutMarker(loc);
                                        currentFamily.getMember(
                                                currentInd.getMomID()).zeroOutMarker(loc);
                                        currentFamily.getMember(
                                                currentInd.getDadID()).zeroOutMarker(loc);
                                    }
                                }
                            }
                        } else {

                            // don't check if parents are missing any data
                            if (!(momAllele1 == 0 || momAllele2 == 0 || dadAllele1 == 0 || dadAllele2 == 0)) {
                                // mom hom
                                if (momAllele1 == momAllele2) {
                                    // both parents hom
                                    if (dadAllele1 == dadAllele2) {
                                        // both parents hom same allele
                                        if (momAllele1 == dadAllele1) {
                                            // kid must be hom same allele
                                            if (allele1 != momAllele1 || allele2 != momAllele1) {
                                                mendErrNum++;
                                                MendelError mend = new MendelError(
                                                        currentInd.getFamilyID(),
                                                        currentInd.getIndividualID());
                                                mendels.add(mend);
                                                currentInd.zeroOutMarker(loc);
                                                currentFamily.getMember(
                                                        currentInd.getMomID()).zeroOutMarker(loc);
                                                currentFamily.getMember(
                                                        currentInd.getDadID()).zeroOutMarker(loc);
                                            }
                                        // parents hom diff allele
                                        } else {
                                            // kid must be het
                                            if (allele1 == allele2) {
                                                mendErrNum++;
                                                MendelError mend = new MendelError(
                                                        currentInd.getFamilyID(),
                                                        currentInd.getIndividualID());
                                                mendels.add(mend);
                                                currentInd.zeroOutMarker(loc);
                                                currentFamily.getMember(
                                                        currentInd.getMomID()).zeroOutMarker(loc);
                                                currentFamily.getMember(
                                                        currentInd.getDadID()).zeroOutMarker(loc);
                                            }
                                        }
                                    // mom hom dad het
                                    } else {
                                        // kid can't be hom for non-momallele
                                        if (allele1 != momAllele1 && allele2 != momAllele1) {
                                            mendErrNum++;
                                            MendelError mend = new MendelError(
                                                    currentInd.getFamilyID(),
                                                    currentInd.getIndividualID());
                                            mendels.add(mend);
                                            currentInd.zeroOutMarker(loc);
                                            currentFamily.getMember(
                                                    currentInd.getMomID()).zeroOutMarker(loc);
                                            currentFamily.getMember(
                                                    currentInd.getDadID()).zeroOutMarker(loc);
                                        }
                                    }
                                // mom het
                                } else {
                                    // dad hom
                                    if (dadAllele1 == dadAllele2) {
                                        // kid can't be hom for non-dadallele
                                        if (allele1 != dadAllele1 && allele2 != dadAllele1) {
                                            mendErrNum++;
                                            MendelError mend = new MendelError(
                                                    currentInd.getFamilyID(),
                                                    currentInd.getIndividualID());
                                            mendels.add(mend);
                                            currentInd.zeroOutMarker(loc);
                                            currentFamily.getMember(
                                                    currentInd.getMomID()).zeroOutMarker(loc);
                                            currentFamily.getMember(
                                                    currentInd.getDadID()).zeroOutMarker(loc);
                                        }
                                    }
                                // both parents het no mend err poss
                                }
                            }
                        }
                    }
                // end mendel check
                }
            }

            indList = currentFamily.getMemberList();
            // loop through each individual in the current Family
            while (indList.hasMoreElements()) {
                currentInd = currentFamily.getMember((String) indList.nextElement());
                if (currentInd.getZeroed(loc)) {
                    allele1 = 0;
                    allele2 = 0;
                } else {
                    allele1 = currentInd.getAllele(loc, 0);
                    allele2 = currentInd.getAllele(loc, 1);
                }

                String familyID = currentInd.getFamilyID();

                // no allele data missing
                if (allele1 > 0 && allele2 > 0) {
                    // indiv has no parents -- i.e. is a founder
                    if (!currentFamily.hasAncestor(currentInd.getIndividualID())) {
                        // set founderGenoCount
                        if (founderGenoCount.containsKey(familyID)) {
                            int value = ((Integer) founderGenoCount.get(familyID)).intValue() + 1;
                            founderGenoCount.put(familyID, new Integer(value));
                        } else {
                            founderGenoCount.put(familyID, new Integer(1));
                        }

                        if (allele1 != 9) { // value of 9 means an 'h' allele for haps files...

                            count[allele1]++;
                        }
                        if (!Chromosome.getDataChrom().equals("chrx") || currentInd.getGender() != 1) {
                            if (allele1 != allele2 || allele1 == 9 || allele2 == 9) {
                                founderHetCount++;
                            } else {
                                founderHomCount[allele1]++;
                            }
                            if (allele2 != 9) {
                                count[allele2]++;
                            }
                        }
                    } else {
                        if (kidgeno.containsKey(familyID)) {
                            int value = ((Integer) kidgeno.get(familyID)).intValue() + 1;
                            kidgeno.put(familyID, new Integer(value));
                        } else {
                            kidgeno.put(familyID, new Integer(1));
                        }
                    }

                    if (!Chromosome.getDataChrom().equals("chrx") || currentInd.getGender() != 1) {
                        if (allele1 == allele2) {
                            hom++;
                        } else {
                            het++;
                        }
                    } else {
                        haploid++;
                    }
                } // missing data
                else {
                    missing++;
                }
            }
            currentFamily.setMendErrs(mendErrNum);
        }
        double obsHET = getObsHET(het, hom);
        double freqStuff[] = null;
        try {
            freqStuff = getFreqStuff(count);
        } catch (MDRPedFileException pfe) {
            throw new MDRPedFileException("More than " + numAllele + " alleles at marker " + (loc + 1));
        }
        double preHET = freqStuff[0];
        double maf = freqStuff[1];
        String minorAllele;
        if (freqStuff[2] == 1) {
            minorAllele = "A";
        } else if (freqStuff[2] == 2) {
            minorAllele = "C";
        } else if (freqStuff[2] == 3) {
            minorAllele = "G";
        } else {
            minorAllele = "T";
        }
        // HW p value
        double pvalue = getPValue(founderHomCount, founderHetCount);

        // geno percent
        double genopct = getGenoPercent(het, hom, haploid, missing);

        // num of families with a fully genotyped trio
        // int famTrio =0;
        int famTrio = getNumOfFamTrio(pedFile.getFamList(), founderGenoCount,
                kidgeno);

        // rating
        int rating = this.getRating(genopct, pvalue, mendErrNum, maf);

        result.setObsHet(obsHET);
        result.setPredHet(preHET);
        result.setMAF(maf);
        result.setMinorAllele(minorAllele);
        result.setHWpvalue(pvalue);
        result.setGenoPercent(genopct);
        result.setFamTrioNum(famTrio);
        result.setMendErrNum(mendErrNum);
        result.setRating(rating);
        result.setMendelErrors(mendels);
        return result;
    }

    /**
     * Gets observed heterozygosity
     */
    private double getObsHET(int het, int hom) {
        double obsHET;
        if (het + hom == 0) {
            obsHET = 0;
        } else {
            obsHET = het / (het + hom + 0.0);
        }
        return obsHET;
    }

    private double[] getFreqStuff(int[] count) throws MDRPedFileException {
        double[] freqStuff = new double[3];
        int sumsq = 0, sum = 0, num = 0, mincount = -1;
        int numberOfAlleles = 0;
        for (int i = 0; i < count.length; i++) {
            if (count[i] != 0) {
                numberOfAlleles++;
                num = count[i];
                sumsq += num * num;
                sum += num;
                if (mincount < 0 || mincount > num) {
                    mincount = num;
                    freqStuff[2] = i;
                }
            }
        }
        /*
         * if (numberOfAlleles > numAllele){ throw new MDRPedFileException("More than two alleles!"); }
         */
        if (sum == 0) {
            freqStuff[0] = 0;
            freqStuff[1] = 0;
        } else {
            freqStuff[0] = 1.0 - (sumsq / ((sum * sum) + 0.0));
            if (mincount / (sum + 0.0) == 1) {
                freqStuff[1] = 0.0;
            } else {
                freqStuff[1] = mincount / (sum + 0.0);
            }
        }
        return freqStuff;
    }

    private double getPValue(int[] parentHom, int parentHet)
            throws MDRPedFileException {
        // ie: 11 13 31 33 -> homA =1 homB = 1 parentHet=2
        int homA = 0, homB = 0;
        double pvalue = 0;
        for (int i = 0; i < parentHom.length; i++) {
            if (parentHom[i] != 0) {
                if (homA > 0) {
                    homB = parentHom[i];
                } else {
                    homA = parentHom[i];
                }
            }
        }
        // caculate p value from homA, parentHet and homB
        if (homA + parentHet + homB <= 0) {
            pvalue = 0;
        } else {
            pvalue = hwCalculate(homA, parentHet, homB);
        }
        return pvalue;
    }

    private double hwCalculate(int obsAA, int obsAB, int obsBB)
            throws MDRPedFileException {
        // Calculates exact two-sided hardy-weinberg p-value. Parameters
        // are number of genotypes, number of rare alleles observed and
        // number of heterozygotes observed.
        //
        // (c) 2003 Jan Wigginton, Goncalo Abecasis
        int diplotypes = obsAA + obsAB + obsBB;
        int rare = (obsAA * 2) + obsAB;
        int hets = obsAB;

        // make sure "rare" allele is really the rare allele
        if (rare > diplotypes) {
            rare = 2 * diplotypes - rare;
        }

        // make sure numbers aren't screwy
        if (hets > rare) {
            throw new MDRPedFileException("HW test: " + hets + "heterozygotes but only " + rare + "rare alleles.");
        }
        double[] tailProbs = new double[rare + 1];
        for (int z = 0; z < tailProbs.length; z++) {
            tailProbs[z] = 0;
        }

        // start at midpoint
        int mid = rare * (2 * diplotypes - rare) / (2 * diplotypes);

        // check to ensure that midpoint and rare alleles have same parity
        if (((rare & 1) ^ (mid & 1)) != 0) {
            mid++;
        }
        int het = mid;
        int hom_r = (rare - mid) / 2;
        int hom_c = diplotypes - het - hom_r;

        // Calculate probability for each possible observed heterozygote
        // count up to a scaling constant, to avoid underflow and overflow
        tailProbs[mid] = 1.0;
        double sum = tailProbs[mid];
        for (het = mid; het > 1; het -= 2) {
            tailProbs[het - 2] = (tailProbs[het] * het * (het - 1.0)) / (4.0 * (hom_r + 1.0) * (hom_c + 1.0));
            sum += tailProbs[het - 2];
            // 2 fewer hets for next iteration -> add one rare and one common homozygote
            hom_r++;
            hom_c++;
        }

        het = mid;
        hom_r = (rare - mid) / 2;
        hom_c = diplotypes - het - hom_r;
        for (het = mid; het <= rare - 2; het += 2) {
            tailProbs[het + 2] = (tailProbs[het] * 4.0 * hom_r * hom_c) / ((het + 2.0) * (het + 1.0));
            sum += tailProbs[het + 2];
            // 2 more hets for next iteration -> subtract one rare and one common homozygote
            hom_r--;
            hom_c--;
        }

        for (int z = 0; z < tailProbs.length; z++) {
            tailProbs[z] /= sum;
        }

        double top = tailProbs[hets];
        for (int i = hets + 1; i <= rare; i++) {
            top += tailProbs[i];
        }
        double otherSide = tailProbs[hets];
        for (int i = hets - 1; i >= 0; i--) {
            otherSide += tailProbs[i];
        }

        if (top > 0.5 && otherSide > 0.5) {
            return 1.0;
        } else {
            if (top < otherSide) {
                return top * 2;
            } else {
                return otherSide * 2;
            }
        }
    }

    private double getGenoPercent(int het, int hom, int haploid, int missing) {
        if (het + hom + haploid + missing == 0) {
            return 0;
        } else {
            return 100.0 * (het + hom + haploid) / (het + hom + haploid + missing);
        }
    }

    private int getNumOfFamTrio(Enumeration famList, Hashtable parentgeno,
            Hashtable kidgeno) {
        // this is buggy. it doesn't do what we want with larger families.
        // it's hard to even define what we want this to represent. oh well.
        int tdtfams = 0;
        while (famList.hasMoreElements()) {
            int parentGeno = 0, kidsGeno = 0;
            String key = (String) famList.nextElement();
            Integer pGeno = (Integer) parentgeno.get(key);
            Integer kGeno = (Integer) kidgeno.get(key);
            if (pGeno != null) {
                parentGeno = pGeno.intValue();
            }
            if (kGeno != null) {
                kidsGeno = kGeno.intValue();
            // basically we want the smaller of either (a) half the number of genotyped parents
            // or (b) the number of genotyped offspring.
            }
            if (parentGeno >= 2 && kidsGeno >= 1) {
                if (parentGeno / 2 > kidsGeno) {
                    tdtfams += kidsGeno;
                } else {
                    tdtfams += parentGeno / 2;
                }
            }
        }
        return tdtfams;
    }

    private int getRating(double genopct, double pval, int menderr, double maf) {
        int rating = 0;
        if (genopct < failedGenoCut) {
            rating -= 2;
        }
        if (pval < hwCut) {
            rating -= 4;
        }
        if (menderr > numMendErrCut) {
            rating -= 8;
        }
        if (maf < mafCut) {
            rating -= 16;
        }
        if (rating == 0) {
            rating = 1;
        }
        return rating;
    }
}
