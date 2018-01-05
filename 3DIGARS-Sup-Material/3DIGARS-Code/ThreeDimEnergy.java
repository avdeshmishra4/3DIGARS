import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;

public class ThreeDimEnergy {

	// alpha and beta's values obtain from 3DIGARS Genetic Algorithm
	// alpha-HP = 1.3802541
	// alpha-HH = 1.6832844
	// alpha-PP = 1.9315737
	// beta-HP = 1.4921875
	// beta-HH = 0.55859375
	// beta-PP = 0.265625
	static float[] alphas = { 1.3802541f, 1.6832844f, 1.9315737f };
	static float[] betas = { 1.4921875f, 0.55859375f, 0.265625f };

	// setup needed for generating and loading frequency table and calculating
	// energy
	static int maxRow = 14028;
	static int maxCol = 30;
	static double small_num = 0.000001;
	static int atomType = 167;

	// frequency table setup

	static float[][] alphaCarbonArray = new float[70000][3];
	static int[] resAtomPairArrayID = new int[70000];
	static String[] pairs = new String[maxRow];
	static int cAlphaCount = -1;
	static int[] resNum = new int[70000];
	static int[] Hphob_Hphil_Array = new int[70000];

	// frequency table setup ends

	static String[][] Hphob_Hphil_Mapper = { { "ALA", "hphob" },
			{ "ARG", "hphil" }, { "ASN", "hphil" }, { "ASP", "hphil" },
			{ "CYS", "hphil" }, { "GLN", "hphil" }, { "GLU", "hphil" },
			{ "GLY", "hphob" }, { "HIS", "hphil" }, { "ILE", "hphob" },
			{ "LEU", "hphob" }, { "LYS", "hphil" }, { "MET", "hphob" },
			{ "PHE", "hphob" }, { "PRO", "hphil" }, { "SER", "hphil" },
			{ "THR", "hphil" }, { "TRP", "hphil" }, { "TYR", "hphob" },
			{ "VAL", "hphob" } }; // this categorization is based on Dr. Hoque's
									// paper

	/*
	 * Possible residue-atom pair
	 */
	static String resAtomPair[] = { "ALA N", "ALA CA", "ALA C", "ALA O",
			"ALA CB", "CYS N", "CYS CA", "CYS C", "CYS O", "CYS CB", "CYS SG",
			"ASP N", "ASP CA", "ASP C", "ASP O", "ASP CB", "ASP CG", "ASP OD1",
			"ASP OD2", "GLU N", "GLU CA", "GLU C", "GLU O", "GLU CB", "GLU CG",
			"GLU CD", "GLU OE1", "GLU OE2", "PHE N", "PHE CA", "PHE C",
			"PHE O", "PHE CB", "PHE CG", "PHE CD1", "PHE CD2", "PHE CE1",
			"PHE CE2", "PHE CZ", "GLY N", "GLY CA", "GLY C", "GLY O", "HIS N",
			"HIS CA", "HIS C", "HIS O", "HIS CB", "HIS CG", "HIS ND1",
			"HIS CD2", "HIS CE1", "HIS NE2", "ILE N", "ILE CA", "ILE C",
			"ILE O", "ILE CB", "ILE CG1", "ILE CG2", "ILE CD1", "LYS N",
			"LYS CA", "LYS C", "LYS O", "LYS CB", "LYS CG", "LYS CD", "LYS CE",
			"LYS NZ", "LEU N", "LEU CA", "LEU C", "LEU O", "LEU CB", "LEU CG",
			"LEU CD1", "LEU CD2", "MET N", "MET CA", "MET C", "MET O",
			"MET CB", "MET CG", "MET SD", "MET CE", "ASN N", "ASN CA", "ASN C",
			"ASN O", "ASN CB", "ASN CG", "ASN OD1", "ASN ND2", "PRO N",
			"PRO CA", "PRO C", "PRO O", "PRO CB", "PRO CG", "PRO CD", "GLN N",
			"GLN CA", "GLN C", "GLN O", "GLN CB", "GLN CG", "GLN CD",
			"GLN OE1", "GLN NE2", "ARG N", "ARG CA", "ARG C", "ARG O",
			"ARG CB", "ARG CG", "ARG CD", "ARG NE", "ARG CZ", "ARG NH1",
			"ARG NH2", "SER N", "SER CA", "SER C", "SER O", "SER CB", "SER OG",
			"THR N", "THR CA", "THR C", "THR O", "THR CB", "THR OG1",
			"THR CG2", "VAL N", "VAL CA", "VAL C", "VAL O", "VAL CB",
			"VAL CG1", "VAL CG2", "TRP N", "TRP CA", "TRP C", "TRP O",
			"TRP CB", "TRP CG", "TRP CD1", "TRP CD2", "TRP NE1", "TRP CE2",
			"TRP CE3", "TRP CZ2", "TRP CZ3", "TRP CH2", "TYR N", "TYR CA",
			"TYR C", "TYR O", "TYR CB", "TYR CG", "TYR CD1", "TYR CD2",
			"TYR CE1", "TYR CE2", "TYR CZ", "TYR OH" };

	static double[][][] hphob_hphob_table = new double[atomType][atomType][maxCol];
	static double[][][] hphil_hphil_table = new double[atomType][atomType][maxCol];
	static double[][][] hphob_hphil_table = new double[atomType][atomType][maxCol];

	/*
	 * holds frequency counts
	 */

	static double[][] freqTbl_hphob_hphob = new double[maxRow][maxCol];
	static double[][] freqTbl_hphil_hphil = new double[maxRow][maxCol];
	static double[][] freqTbl_hphob_hphil = new double[maxRow][maxCol];

	/*
	 * holds energy scores
	 */

	static double[][][] probTblHphob = new double[atomType][atomType][maxCol];
	static double[][][] probTblHphil = new double[atomType][atomType][maxCol];
	static double[][][] probTblHphob_Hphil = new double[atomType][atomType][maxCol];

	String[] pairID = new String[maxRow]; // this will carry the string of atom
											// pair --- here [8] is to hold the
											// array of at most (167-167) seven
											// characters

	// this variable stores the energy value for given structure
	static double predictedEnergy = 0.;

	public static void main(String[] args) {

		try {

			File currentDir = new File(new File(".").getAbsolutePath());
			String curr_dir_path = currentDir.getCanonicalPath();
			
			if(args.length == 0){
				
				System.out.println("Please provide the full path and file name of the structure for which you want to calculate energy");
				System.exit(0);
			}
			
			String inputFile = args[0];

			/*
			 * uncomment the following line if you want to generate the energy score libraries again
			 */
			
			// generateFreqAndProbTable(alphas[0], alphas[1], alphas[2],curr_dir_path);

			readProbTableInMemory(curr_dir_path);			// read energy score library in memory

			predictedEnergy = 0;

			findCAlphaAndLoadInMemory(inputFile);			// parses all the ATOMS line from pdb structure and puts in memory

			calcEucDistAndProtability(betas[0], betas[1], betas[2]);	// calculates the energy score for given structure

			float ThreeDIGARSEnergy = (float) (predictedEnergy / 100);	// scale the final output by dividing it with 100

			System.out.println(inputFile + "\t" + ThreeDIGARSEnergy);	// print the final result

		} catch (Exception e) {

			e.printStackTrace();

		}

	}

	public static void generateFreqAndProbTable(float a1, float a2, float a3,
			String curr_dir_path) throws IOException {

		String pdbDir = curr_dir_path + "/3DIGARSTrainingDataPDB/";

		String probTblStr_hphob_hphob = curr_dir_path + "/probTableHphob.txt";

		String probTblStr_hphil_hphil = curr_dir_path + "/probTableHphil.txt";

		String probTblStr_hphob_hphil = curr_dir_path
				+ "/probTableHphob_Hphil.txt";

		File probTbl_file_hphob_hphob = new File(probTblStr_hphob_hphob);
		File probTbl_file_hphil_hphil = new File(probTblStr_hphil_hphil);
		File probTbl_file_hphob_hphil = new File(probTblStr_hphob_hphil);

		BufferedWriter prob_hphob_hphob_writer = BufferReaderAndWriter
				.getWriter(probTbl_file_hphob_hphob);
		BufferedWriter prob_hphil_hphil_writer = BufferReaderAndWriter
				.getWriter(probTbl_file_hphil_hphil);
		BufferedWriter prob_hphob_hphil_writer = BufferReaderAndWriter
				.getWriter(probTbl_file_hphob_hphil);

		String fastaInput = curr_dir_path + "/3DIGARSTrainDataFasta.txt";

		freqTblFirstColFiller();

		File fastaFile = new File(fastaInput);
		BufferedReader fastaReader = BufferReaderAndWriter.getReader(fastaFile);
		String fastaLine = "";
		while ((fastaLine = fastaReader.readLine()) != null) {

			if (fastaLine.startsWith(">")) {
				cAlphaCount = -1;

				// System.out.println("FastaID => " + fastaLine);

				String fileName = fastaLine.substring(1, 5).toLowerCase();

				String fileNameWithDir = pdbDir + fileName + ".ent";

				findCAlphaAndLoadInMemory(fileNameWithDir);

				calcEucDistAndUpdateFreqTable();

			}

		}

		fastaReader.close();

		int counter1 = -1;
		for (int i = 0; i < atomType; i++) {

			for (int j = i; j < atomType; j++) {
				counter1++;
				for (int k = 0; k < maxCol; k++) {

					freqTbl_hphob_hphob[counter1][k] = hphob_hphob_table[i][j][k];

				}

			}

		}

		// convert 3 dim array into 2 dim array for hphil_hphil

		int counter2 = -1;
		for (int i = 0; i < atomType; i++) {

			for (int j = i; j < atomType; j++) {
				counter2++;
				for (int k = 0; k < maxCol; k++) {

					freqTbl_hphil_hphil[counter2][k] = hphil_hphil_table[i][j][k];

				}

			}

		}

		// convert 3 dim array into 2 dim array for hphob_hphil or hphil_hphob
		// (symmetric)

		int counter3 = -1;
		for (int i = 0; i < atomType; i++) {

			for (int j = i; j < atomType; j++) {
				counter3++;
				for (int k = 0; k < maxCol; k++) {

					freqTbl_hphob_hphil[counter3][k] = hphob_hphil_table[i][j][k];

				}

			}

		}

		// // replace zero with a small value for hphob_hphob

		for (int x = 0; x < maxRow; x++) {

			for (int j = 0; j < maxCol; j++) {

				if (freqTbl_hphob_hphob[x][j] == 0.0) {

					freqTbl_hphob_hphob[x][j] = small_num;

				}

			}

		}

		// replace zero with a small value for hphil_hphil

		for (int x = 0; x < maxRow; x++) {

			for (int j = 0; j < maxCol; j++) {

				if (freqTbl_hphil_hphil[x][j] == 0.0) {

					freqTbl_hphil_hphil[x][j] = small_num;

				}

			}

		}

		// replace zero with a small value for hphob_hphil or hphil_hphob

		for (int x = 0; x < maxRow; x++) {

			for (int j = 0; j < maxCol; j++) {

				if (freqTbl_hphob_hphil[x][j] == 0.0) {

					freqTbl_hphob_hphil[x][j] = small_num;

				}

			}

		}

		// generate energy score matrix and print in a file for hphob
		double r_phob = 0.;
		double r_cut_phob = 15.;
		double delta_r_phob = 0.5;
		double delta_r_cut_phob = 0.5;
		// double alpha_phob = 1.;
		double nr_phob = 0.; // numerator
		double dr_phob = 0.; // denomenator
		double element_phob = 0.;

		for (int t = 0; t < maxRow; t++) {
			r_phob = 0.;
			prob_hphob_hphob_writer.write(pairs[t]);
			prob_hphob_hphob_writer.write(",");

			for (int j = 0; j < maxCol; j++) {

				if (j < maxCol - 1) {
					r_phob = r_phob + 0.5;

					nr_phob = freqTbl_hphob_hphob[t][j];

					double lastElementOfAllTable = freqTbl_hphob_hphob[t][maxCol - 1]
							+ freqTbl_hphil_hphil[t][maxCol - 1]
							+ freqTbl_hphob_hphil[t][maxCol - 1];
					dr_phob = java.lang.Math.pow((r_phob / r_cut_phob), a1)
							* (delta_r_phob / delta_r_cut_phob)
							* (lastElementOfAllTable);
					element_phob = -java.lang.Math.log(nr_phob / dr_phob);

					prob_hphob_hphob_writer
							.write(Double.toString(element_phob));
					prob_hphob_hphob_writer.write(",");
					prob_hphob_hphob_writer.flush();

				} else if (j == maxCol - 1) {
					r_phob = r_phob + 0.5;

					nr_phob = freqTbl_hphob_hphob[t][j];
					double lastElementOfAllTable = freqTbl_hphob_hphob[t][maxCol - 1]
							+ freqTbl_hphil_hphil[t][maxCol - 1]
							+ freqTbl_hphob_hphil[t][maxCol - 1];
					dr_phob = java.lang.Math.pow((r_phob / r_cut_phob), a1)
							* (delta_r_phob / delta_r_cut_phob)
							* (lastElementOfAllTable);
					element_phob = -java.lang.Math.log(nr_phob / dr_phob);

					prob_hphob_hphob_writer
							.write(Double.toString(element_phob));
					prob_hphob_hphob_writer.newLine();
					prob_hphob_hphob_writer.flush();

				}

			}

		}

		prob_hphob_hphob_writer.close();

		// generate energy score matrix and print in a file for hphil
		double r_hphil = 0.;
		double r_cut_hphil = 15.;
		double delta_r_hphil = 0.5;
		double delta_r_cut_hphil = 0.5;
		// double alpha_hphil = 2.;
		double nr_hphil = 0.; // numerator
		double dr_hphil = 0.; // denomenator
		double element_hphil = 0.;

		for (int t = 0; t < maxRow; t++) {
			r_hphil = 0.;
			prob_hphil_hphil_writer.write(pairs[t]);
			prob_hphil_hphil_writer.write(",");

			for (int j = 0; j < maxCol; j++) {

				if (j < maxCol - 1) {
					r_hphil = r_hphil + 0.5;
					// nr = frequencyTable[t][j] / rowSum[t];
					nr_hphil = freqTbl_hphil_hphil[t][j];
					// dr = colSum[j] / allColSum;
					double lastElementOfAllTable = freqTbl_hphob_hphob[t][maxCol - 1]
							+ freqTbl_hphil_hphil[t][maxCol - 1]
							+ freqTbl_hphob_hphil[t][maxCol - 1];
					dr_hphil = java.lang.Math.pow((r_hphil / r_cut_hphil), a2)
							* (delta_r_hphil / delta_r_cut_hphil)
							* (lastElementOfAllTable);
					element_hphil = -java.lang.Math.log(nr_hphil / dr_hphil);

					prob_hphil_hphil_writer.write(Double
							.toString(element_hphil));
					prob_hphil_hphil_writer.write(",");
					prob_hphil_hphil_writer.flush();

				} else if (j == maxCol - 1) {
					r_hphil = r_hphil + 0.5;
					// nr = frequencyTable[t][j] / rowSum[t];
					nr_hphil = freqTbl_hphil_hphil[t][j];
					// dr = colSum[j] / allColSum;
					double lastElementOfAllTable = freqTbl_hphob_hphob[t][maxCol - 1]
							+ freqTbl_hphil_hphil[t][maxCol - 1]
							+ freqTbl_hphob_hphil[t][maxCol - 1];
					dr_hphil = java.lang.Math.pow((r_hphil / r_cut_hphil), a2)
							* (delta_r_hphil / delta_r_cut_hphil)
							* (lastElementOfAllTable);
					element_hphil = -java.lang.Math.log(nr_hphil / dr_hphil);

					prob_hphil_hphil_writer.write(Double
							.toString(element_hphil));
					prob_hphil_hphil_writer.newLine();
					prob_hphil_hphil_writer.flush();

				}

			}

		}

		prob_hphil_hphil_writer.close();

		// generate energy score matrix and print in a file for hphob_hphil
		double r_hphob_hphil = 0.;
		double r_cut_hphob_hphil = 15.;
		double delta_r_hphob_hphil = 0.5;
		double delta_r_cut_hphob_hphil = 0.5;
		// double alpha_hphob_hphil = 1.57;
		double nr_hphob_hphil = 0.; // numerator
		double dr_hphob_hphil = 0.; // denomenator
		double element_hphob_hphil = 0.;

		for (int t = 0; t < maxRow; t++) {
			r_hphob_hphil = 0.;

			prob_hphob_hphil_writer.write(pairs[t]);
			prob_hphob_hphil_writer.write(",");

			for (int j = 0; j < maxCol; j++) {

				if (j < maxCol - 1) {
					r_hphob_hphil = r_hphob_hphil + 0.5;
					// nr = frequencyTable[t][j] / rowSum[t];
					nr_hphob_hphil = freqTbl_hphob_hphil[t][j];
					// dr = colSum[j] / allColSum;
					double lastElementOfAllTable = freqTbl_hphob_hphob[t][maxCol - 1]
							+ freqTbl_hphil_hphil[t][maxCol - 1]
							+ freqTbl_hphob_hphil[t][maxCol - 1];
					dr_hphob_hphil = java.lang.Math.pow(
							(r_hphob_hphil / r_cut_hphob_hphil), a3)
							* (delta_r_hphob_hphil / delta_r_cut_hphob_hphil)
							* (lastElementOfAllTable);
					element_hphob_hphil = -java.lang.Math.log(nr_hphob_hphil
							/ dr_hphob_hphil);

					prob_hphob_hphil_writer.write(Double
							.toString(element_hphob_hphil));
					prob_hphob_hphil_writer.write(",");
					prob_hphob_hphil_writer.flush();

				} else if (j == maxCol - 1) {
					r_hphob_hphil = r_hphob_hphil + 0.5;
					// nr = frequencyTable[t][j] / rowSum[t];
					nr_hphob_hphil = freqTbl_hphob_hphil[t][j];
					// dr = colSum[j] / allColSum;
					double lastElementOfAllTable = freqTbl_hphob_hphob[t][maxCol - 1]
							+ freqTbl_hphil_hphil[t][maxCol - 1]
							+ freqTbl_hphob_hphil[t][maxCol - 1];
					dr_hphob_hphil = java.lang.Math.pow(
							(r_hphob_hphil / r_cut_hphob_hphil), a3)
							* (delta_r_hphob_hphil / delta_r_cut_hphob_hphil)
							* (lastElementOfAllTable);
					element_hphob_hphil = -java.lang.Math.log(nr_hphob_hphil
							/ dr_hphob_hphil);

					prob_hphob_hphil_writer.write(Double
							.toString(element_hphob_hphil));
					prob_hphob_hphil_writer.newLine();
					prob_hphob_hphil_writer.flush();

				}

			}

		}

		prob_hphob_hphil_writer.close();

		// System.out.println("Process Completed");

	}

	static void calcEucDistAndUpdateFreqTable() {

		for (int x = 0; x < cAlphaCount + 1; x++) {

			for (int y = x + 1; y < cAlphaCount + 1; y++) {

				if (resNum[x] == resNum[y])
					continue;

				int col = 0;

				double eucDis = 0.;

				for (int i = 0; i < 3; i++) {

					double sqr = alphaCarbonArray[x][i]
							- alphaCarbonArray[y][i];

					eucDis += sqr * sqr;

				}

				eucDis = java.lang.Math.sqrt(eucDis);
				col = (int) (eucDis * 2);

				if (col >= maxCol)
					continue;

				if (Hphob_Hphil_Array[x] == 1 && Hphob_Hphil_Array[y] == 1) {
					hphob_hphob_table[resAtomPairArrayID[x]][resAtomPairArrayID[y]][col] += 1;
					hphob_hphob_table[resAtomPairArrayID[y]][resAtomPairArrayID[x]][col] = hphob_hphob_table[resAtomPairArrayID[x]][resAtomPairArrayID[y]][col];
				} else if (Hphob_Hphil_Array[x] == 0
						&& Hphob_Hphil_Array[y] == 0) {
					hphil_hphil_table[resAtomPairArrayID[x]][resAtomPairArrayID[y]][col] += 1;
					hphil_hphil_table[resAtomPairArrayID[y]][resAtomPairArrayID[x]][col] = hphil_hphil_table[resAtomPairArrayID[x]][resAtomPairArrayID[y]][col];

				} else if ((Hphob_Hphil_Array[x] == 1 && Hphob_Hphil_Array[y] == 0)
						|| (Hphob_Hphil_Array[x] == 0 && Hphob_Hphil_Array[y] == 1)) {

					hphob_hphil_table[resAtomPairArrayID[x]][resAtomPairArrayID[y]][col] += 1;
					hphob_hphil_table[resAtomPairArrayID[y]][resAtomPairArrayID[x]][col] = hphob_hphil_table[resAtomPairArrayID[x]][resAtomPairArrayID[y]][col];

				}

			}

		}

	}

	static void findCAlphaAndLoadInMemory(String fileName) throws IOException {

		File file = new File(fileName);
	
	System.out.println(fileName);
		BufferedReader fileReader = BufferReaderAndWriter.getReader(file);

		String atomName = null;
		String resName = null;

		String resCheck = "UNK";
		String distRes = null;
		int resID = -1;
		String ignore = null;

		String pdbFileLine = "";

		while ((pdbFileLine = fileReader.readLine()) != null) {

			if (!pdbFileLine.startsWith("ATOM")) {
				continue;
			}

			atomName = pdbFileLine.substring(13, 16).trim();
			if (atomName.startsWith("H"))
				continue;

			resName = pdbFileLine.substring(17, 20).trim();

			ignore = resName + " " + atomName; // to check if residue name and
												// atom name pair is a valid
												// pair based on the resAtomPair
												// decleared above

			int found = 0;
			for (int i = 0; i < atomType; i++) {

				if (ignore.equals(resAtomPair[i])) {
					cAlphaCount++;
					found = 1;
					resAtomPairArrayID[cAlphaCount] = i;
					break;

				}

			}

			if (found == 0) {

				continue;
			}

			distRes = pdbFileLine.substring(17, 26).trim();

			if (!resCheck.equals(distRes)) {

				resCheck = distRes;

				resID++;

			}

			for (int i = 0; i < 20; i++) {

				if (Hphob_Hphil_Mapper[i][0].equals(resName)) {

					if (Hphob_Hphil_Mapper[i][1].equals("hphob")) {

						Hphob_Hphil_Array[cAlphaCount] = 1;

					} else if (Hphob_Hphil_Mapper[i][1].equals("hphil")) {

						Hphob_Hphil_Array[cAlphaCount] = 0;

					}

					break;

				}

			}

			// printf("%d\n", Hphob_Hphil_Array[cAlphaCount]);

			resNum[cAlphaCount] = resID;

			String xCor = pdbFileLine.substring(30, 38).trim();
			String yCor = pdbFileLine.substring(38, 46).trim();
			String zCor = pdbFileLine.substring(46, 54).trim();

			alphaCarbonArray[cAlphaCount][0] = Float.parseFloat(xCor);
			alphaCarbonArray[cAlphaCount][1] = Float.parseFloat(yCor);
			alphaCarbonArray[cAlphaCount][2] = Float.parseFloat(zCor);

		}

		fileReader.close();
		// printf("done reading file %s\n", fileName);
	}

	static void freqTblFirstColFiller() {

		int count = -1;

		for (int i = 0; i < atomType; i++) {

			for (int j = i; j < atomType; j++) {
				count++;

				String pair = resAtomPair[i] + "," + resAtomPair[j];
				pairs[count] = pair;

			}

		}

	}

	public static void readProbTableInMemory(String curr_dir_path)
			throws Exception {
		String probTableHphob = curr_dir_path + "/probTableHphob.txt";
		System.out.println(probTableHphob);
		String probTableHphil = curr_dir_path + "/probTableHphil.txt";
		String probTableHphob_Hphil = curr_dir_path
				+ "/probTableHphob_Hphil.txt";
		File file_hphob = new File(probTableHphob);
		BufferedReader br_hphob = BufferReaderAndWriter.getReader(file_hphob);
		File file_hphil = new File(probTableHphil);
		BufferedReader br_hphil = BufferReaderAndWriter.getReader(file_hphil);
		File file_hphob_hphil = new File(probTableHphob_Hphil);
		BufferedReader br_hphob_hphil = BufferReaderAndWriter
				.getReader(file_hphob_hphil);

		// load energy score in memory (probTableHphob.txt)
		String resAHphob;
		String resBHphob;

		String line_hphob = "";
		while ((line_hphob = br_hphob.readLine()) != null) {

			if (line_hphob.startsWith("#"))
				continue;

			String tokenHphob[] = line_hphob.split(",");

			int resIDA = -1;
			int resIDB = -1;

			resAHphob = tokenHphob[0];
			resBHphob = tokenHphob[1];

			for (int i = 0; i < atomType; i++) {

				if (resAHphob.equals(resAtomPair[i])) {

					resIDA = i;
				}

				if (resBHphob.equals(resAtomPair[i])) {

					resIDB = i;

				}

				if (resIDA >= atomType || resIDB >= atomType) {
					System.out
							.println("something wrong while loading probability table");
					System.exit(1);
				}

				if (resIDA > -1 && resIDB > -1) {

					break;

				}

			}

			for (int j = 2; j < tokenHphob.length; j++) {

				probTblHphob[resIDA][resIDB][j - 2] = Double
						.parseDouble(tokenHphob[j]);
				probTblHphob[resIDB][resIDA][j - 2] = Double
						.parseDouble(tokenHphob[j]);
			}

		}
		br_hphob.close();

		// load energy score in memory (probTableHphil.txt)

		resAHphob = null;
		resBHphob = null;

		String line_hphil = "";

		while ((line_hphil = br_hphil.readLine()) != null) {

			if (line_hphil.startsWith("#"))
				continue;

			String tokenHphil[] = line_hphil.split(",");

			int resIDA = -1;
			int resIDB = -1;

			resAHphob = tokenHphil[0];
			resBHphob = tokenHphil[1];

			for (int i = 0; i < atomType; i++) {

				if (resAHphob.equals(resAtomPair[i])) {

					resIDA = i;
				}

				if (resBHphob.equals(resAtomPair[i])) {

					resIDB = i;

				}

				if (resIDA >= atomType || resIDB >= atomType) {
					System.out
							.println("something wrong while loading probability table");
					System.exit(1);
				}

				if (resIDA > -1 && resIDB > -1) {

					break;

				}

			}

			for (int j = 2; j < tokenHphil.length; j++) {

				probTblHphil[resIDA][resIDB][j - 2] = Double
						.parseDouble(tokenHphil[j]);
				probTblHphil[resIDB][resIDA][j - 2] = Double
						.parseDouble(tokenHphil[j]);
			}

		}

		br_hphil.close();

		// load energy score in memory (probTableHphob_Hphil.txt)

		resAHphob = null;
		resBHphob = null;

		String line_hphob_hphil = "";

		while ((line_hphob_hphil = br_hphob_hphil.readLine()) != null) {

			if (line_hphob_hphil.startsWith("#"))
				continue;

			String token_Hphob_Hphil[] = line_hphob_hphil.split(",");

			int resIDA = -1;
			int resIDB = -1;

			resAHphob = token_Hphob_Hphil[0];
			resBHphob = token_Hphob_Hphil[1];

			for (int i = 0; i < atomType; i++) {

				if (resAHphob.equals(resAtomPair[i])) {

					resIDA = i;
				}

				if (resBHphob.equals(resAtomPair[i])) {

					resIDB = i;

				}

				if (resIDA >= atomType || resIDB >= atomType) {
					System.out
							.println("something wrong while loading probability table");
					System.exit(1);
				}

				if (resIDA > -1 && resIDB > -1) {

					break;

				}

			}

			for (int j = 2; j < token_Hphob_Hphil.length; j++) {

				probTblHphob_Hphil[resIDA][resIDB][j - 2] = Double
						.parseDouble(token_Hphob_Hphil[j]);
				probTblHphob_Hphil[resIDB][resIDA][j - 2] = Double
						.parseDouble(token_Hphob_Hphil[j]);
			}

		}

		br_hphob_hphil.close();

		// end of reading library file

	}

	public static void calcEucDistAndProtability(float b1, float b2, float b3) {

		predictedEnergy = 0.;
		double pred_hphob = 0.;
		double pred_hphil = 0.;
		double pred_hphob_hphil = 0.;

		for (int x = 0; x < cAlphaCount + 1; x++) {

			for (int y = x + 1; y < cAlphaCount + 1; y++) {

				if (resNum[x] == resNum[y])
					continue;

				int col = 0;

				double eucDis = 0.;

				for (int i = 0; i < 3; i++) {

					double sqr = alphaCarbonArray[x][i]
							- alphaCarbonArray[y][i];

					eucDis += sqr * sqr;

				}

				eucDis = java.lang.Math.sqrt(eucDis);
				col = (int) (eucDis * 2);

				if (col >= maxCol)
					continue;

				if (Hphob_Hphil_Array[x] == 1 && Hphob_Hphil_Array[y] == 1) {

					pred_hphob += probTblHphob[resAtomPairArrayID[x]][resAtomPairArrayID[y]][col];

				} else if (Hphob_Hphil_Array[x] == 0
						&& Hphob_Hphil_Array[y] == 0) {

					pred_hphil += probTblHphil[resAtomPairArrayID[x]][resAtomPairArrayID[y]][col];

				} else if ((Hphob_Hphil_Array[x] == 1 && Hphob_Hphil_Array[y] == 0)
						|| (Hphob_Hphil_Array[x] == 0 && Hphob_Hphil_Array[y] == 1)) {

					pred_hphob_hphil += probTblHphob_Hphil[resAtomPairArrayID[x]][resAtomPairArrayID[y]][col];

				}

			}

		}

		predictedEnergy = (b1 * pred_hphob) + (b2 * pred_hphil)
				+ (b3 * pred_hphob_hphil);

		// printf("predicted value %lf\n", predictedEnergy);

	}

}
