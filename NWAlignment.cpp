/*	This programme alignes DNA/RNA or protein sequences
	provided by the user based on the Needleman-Wunsch algorithm.
	Written by Kyran Wissink on 12/12/2020 */


	/************************************
	*			  Libraries				*
	*************************************/

#include <iostream> // clogs and couts
#include <fstream> // in and out filestream
#include <cstdlib> // exit codes
#include <string> // strings
#include <vector> // vectors
#include <algorithm> // vector and string manipulation
#include <chrono> // duration of code
#include "PAM250.h" // the similarity matrix



/************************************
*			Global variables		*
*************************************/

// Data input and output variables
std::vector <std::string> vecSequences; // input sequences
std::vector <std::string> vecFastaInfo; // info about fasta files
std::string strFileType; // the filetype of the input
bool bProtein = false; // true if aligning protein sequences

std::ofstream ofs("AlignmentResult.txt"); // global filestream access
int iAlignment = 1; // to correctly name alignments in output
int iIdentityCount = 0; // amount of identical amino acids
int iSimilarityCount = 0;// amount of similar amino acids
int iGapCount = 0; // amount of gaps

// Needleman-Wunsch table variables
const double iMatch = 1;
const double iMisMatch = -2;
const double gapOpen = -10;
const double gap = -0.5;

int seq0 = 0; // sequence indices used in the table
int seq1 = 1; // makes multiple sequences technically possible
int x, y; // coordinates on table
int xAxis, yAxis; // size of the axes based on the sequence lengths
std::vector<std::vector <double> > NWtable; // Needleman-Wunsch table
std::vector<std::vector <std::string> > alignmentTracker; // stores routes


// alignment structure
struct alt { // to track alternative routes of alignments
	int x; // that are not yet taken
	int y;
	char route;
};

// alignment variables
int alignpercent = 0; // tracks the percentage of alignment
int index0 = 0; // alignment indices
int index1 = 1;
std::vector <std::string> vecAlignment(2); // The main alignment
char cDirection; // stores the chosen route

//alternative variables
const int maxAltCount = 100000; // maximum number of alts (>500k not recommended)
int altCount = 0; // current number of alts calculated
std::vector <alt> altTracker; // stores alternative routes



/************************************
*			Data reading			*
*************************************/

// checks the validity of sequences and corrects when necessary
std::string validityChecker(std::string& sequence)
{
	// DNA base pairs and amino acids
	std::string strLegalChars = "ATCGU-BDEFGHIJKLMNOPQRSTUVWXYZ*";
	std::string newSequence;

	// protein check
	std::string strNucleotides = "ATGCU-";

	// loop for correcting nucleotide mistakes
	for (int i = 0; i < sequence.size(); ++i) {
		char cSymbol = sequence[i];

		// checks once if it is a protein
		if (bProtein == false && 
			strNucleotides.find(cSymbol) == std::string::npos) {
			bProtein = true;
		}

		// only add the character if it is a valid character
		if (strLegalChars.find(cSymbol) != std::string::npos)
			newSequence += cSymbol;

		// output to log if an illegal character was skipped
		else
			std::clog << "Illegal character skipped: \"" << cSymbol
			<< "\" at position: " << i << '\n';
	}
	return newSequence;
}



// reads the data 
void dataReader(std::string& dataFile) {
	try {

		// open filestream
		std::ifstream ifs(dataFile);
		if (!ifs.is_open())
			throw std::runtime_error("unable to open file.\n"
				"Please check the file path.\n");
		else
			std::clog << "\nReading data from file: " << dataFile << "\n";

		// check filetype by removing the filename
		strFileType = dataFile;
		while (strFileType.front() != '.')
			strFileType.erase(strFileType.begin());

		// fasta files can have multiple extensions
		if (strFileType[1] == 'f') strFileType = ".fasta";

		// txt format (mostly for debugging)
		if (strFileType == ".txt") {
			std::string line, strSequence;
			do {
				ifs >> strSequence;
				strSequence = validityChecker(strSequence);

				// check if the sequence is not empty
				if (!strSequence.empty())
					vecSequences.push_back(strSequence);
				else
					throw std::runtime_error("invalid DNA sequence given.\n");
			} while (std::getline(ifs, line));
		}

		// all FASTA formats
		else if (strFileType == ".fasta") {
			std::string info, strSequence, fullSequence;
			std::getline(ifs, info); // first line is always sequence information
			vecFastaInfo.push_back(info); // store info
			while (ifs.good()) {
				ifs >> strSequence;
				if (strSequence.front() == '>') { // indicates new sequence
					fullSequence = validityChecker(fullSequence);
					vecSequences.push_back(fullSequence); // stores sequence
					fullSequence.clear();
					std::getline(ifs, info);
					strSequence.append(info);
					info = strSequence; // stores next info
					vecFastaInfo.push_back(info);
				}
				else // sequences span multiple lines in FASTA
					fullSequence.append(strSequence);
			}
			fullSequence = validityChecker(fullSequence);
			vecSequences.push_back(fullSequence);
		}

		else {
			throw std::invalid_argument
			("unsupported filetype specified in arguments.\n");
		}

		ifs.close();
		if (ifs.is_open()) {
			throw std::runtime_error("unable to close file.\n");
		}
		else
			std::clog << "Succesfully finished reading: " << dataFile << "\n";
	}

	catch (std::exception& error) {
		std::cerr << "\nError: " << error.what();
		exit(EXIT_FAILURE);
	}
}



/************************************
*		Needleman-Wunsch table		*
*************************************/

// computes the NW table
void tableMaker()
{
	// Initialisation
	std::clog << "Initialising...\n";

	// axes are 1 over the sequence length, since the first row and column  
	// are independent of the sequences and always the gap score
	xAxis = vecSequences[seq0].size() + 1;
	yAxis = vecSequences[seq1].size() + 1;
	NWtable = std::vector<std::vector <double> >
		(xAxis, std::vector<double>(yAxis));
	alignmentTracker = std::vector<std::vector <std::string> >
		(xAxis, std::vector<std::string>(yAxis));
	int tablePercent = 0; // tracks completion of the table

	// sets the first row and column to be proportionate to the gap score.
	for (int x = 0; x < xAxis; ++x) {
		NWtable[x][0] = x * gap;
		alignmentTracker[x][0] = 'l';
	}
	for (int y = 0; y < yAxis; ++y) {
		NWtable[0][y] = y * gap;
		alignmentTracker[x][0] = 'u';
	}

	// initialise algorithm scores
	double uScore, lScore, dScore, mScore;

	// Main algorithm
	for (int x = 1, y = 1; y < yAxis; ++x) {
		if (x == xAxis) { // if the end of the table is reached, move up
			++y;
			x = 0;
			continue;
		}

		// check the diagonal score
		char c1 = vecSequences[seq0][x - 1];
		char c2 = vecSequences[seq1][y - 1];

		// look up the score on the PAM250 matrix for proteins
		if (bProtein == true) {
			int it1 = itAA.find(c1);
			int it2 = itAA.find(c2);
			mScore = PAM250[it1][it2];
		}
		else { // for DNA/RNA, simply score match or mismatch
			if (c1 == c2)
				mScore = iMatch;
			else
				mScore = iMisMatch;
		}

		// gets the highest score from all possible cells to determine score
		if (alignmentTracker[x][y - 1][0] != 'u' && y != yAxis - 1)
			uScore = NWtable[x][y - 1] + gapOpen; // score from upper cell
		else
			uScore = NWtable[x][y - 1] + gap;
		if (alignmentTracker[x - 1][y][0] != 'l' && x != xAxis - 1)
			lScore = NWtable[x - 1][y] + gapOpen; // score from left cell
		else
			lScore = NWtable[x - 1][y] + gap;
		dScore = NWtable[x - 1][y - 1] + mScore; // score from diag cell
		double max[] = { uScore, lScore, dScore };
		double maxScore = *std::max_element(max, max + 3); // highest score
		NWtable[x][y] = maxScore; // input the highest score in the table

		// store equally viable routes for alternative alignments
		std::string routes;
		if (dScore == maxScore)
			routes.push_back('d');
		if (uScore == maxScore)
			routes.push_back('u');
		if (lScore == maxScore)
			routes.push_back('l');

		// add possible routes to the alignment tracker (for alts)
		alignmentTracker[x][y] = routes;

		// calculate current percentage based on the calculated coordinates 
		int tableNewPercent = y / ((static_cast<double> (yAxis)) - 1) * 100;
		if (tableNewPercent != tablePercent) {
			std::cout << "Constructing the Needleman-Wunsch table: "
				<< tableNewPercent << "%\r";
			tablePercent = tableNewPercent;
		}
	}
	std::clog << "\nSuccesfully constructed Needleman-Wunsch table!\n\n";
}



// shows the NWtable on the screen (useful in debugging)
void showTable()
{
	std::cout << "\nNeedleman-Wunsch table\n\t";

	// prints the nucleotides on the top
	for (int i = 0; i < vecSequences[seq0].size(); ++i) {
		std::cout << vecSequences[seq0][i] << "   ";
	}
	std::cout << "\n ";

	// prints the rest of the table
	for (int x = 0, y = 0, spacecount = 4; y < yAxis; ++x, spacecount = 3) {

		// the other sequence should be on the y-axis
		if (x == xAxis) { // if the end of the axis is reached, move up and
			std::cout << '\n' << vecSequences[seq1][y]; // print next 
			++y;						 // character of the other sequence
			x = -1;
			continue;
		}

		// amount of spaces depend on character length of the integer
		if (NWtable[x][y] / 10 != 0)
			spacecount -= 1; // one less space if the value is over 10
		if (NWtable[x][y] < 0)
			spacecount -= 1; // one less space if its minus

		// input spaces and then the value
		for (int i = 0; i < spacecount; ++i)
			std::cout << " ";
		std::cout << NWtable[x][y];
	}
}



/************************************
*		Conserativeness check		*
*************************************/

char ConserveCheck(char& c1, char& c2)
{
	// if there is a gap
	if (c1 == '-' || c2 == '-')
		return ' ';

	// look up the score on the PAM250 matrix for similarity score
	int it1 = itAA.find(c1); // gets row iterator for amino acid
	int it2 = itAA.find(c2); // gets column iterator for amino acid
	if (PAM250[it1][it2] > 0) 
		return ':'; // bigger than 0 = relatively similar
	else
		return '.'; // not very similar
}



/************************************
*			Output stream			*
*************************************/

// checks the identity, similarity and gap scores
void scoreChecker()
{
	for (int i = 0; i < vecAlignment[index0].size(); ++i) {
		char c1 = vecAlignment[index0][i];
		char c2 = vecAlignment[index1][i];

		// if the nucleotides/amino acids are the same
		if (c1 == c2) {
			++iIdentityCount;
			continue;
		}

		// test similarity of amino acids
		if (bProtein == true) {
			char simScore = ConserveCheck(c1, c2);
			if (simScore == ' ')
				++iGapCount;
			else if (simScore == ':')
				++iSimilarityCount;
			// no score if '.' because it is not similar
		}
	}
}



// makes a header for the output file
void outputHeader()
{
	// calculate percentages
	float identitypercent = iIdentityCount /
		static_cast<float>(vecAlignment[index0].size()) * 100;
	float gappercent = iGapCount /
		static_cast<float>(vecAlignment[index0].size()) * 100;
	float similaritypercent = (iSimilarityCount + iIdentityCount) /
		static_cast<float>(vecAlignment[index0].size()) * 100;

	// begin of header
	ofs << "#############################################################\n";
	ofs << "Results of Needleman-Wunsch alignment of:\n";

	// print fasta info if possible, otherwise print sequence
	if (strFileType == ".fasta") {
		for (int i = 0; i < vecFastaInfo.size(); ++i)
			ofs << vecFastaInfo[i] << '\n';
	}
	else {
		for (int i = 0; i < vecSequences.size() - 1; ++i)
			ofs << vecSequences[i] << '\n';
	}

	// print the scoring system
	ofs << "\nScoring system\n"; 
	
	// nucleotides
	if (bProtein == false) {
		ofs << "Match score :\t" << iMatch << '\n';
		ofs << "Mismatch score:\t" << iMisMatch << '\n';
	}
	else { // proteins
		ofs << "Similarity matrix used:\tPAM250" << '\n';
	}

	ofs << "Gap create penalty:\t" << gapOpen << '\n';
	ofs << "Gap extend penalty:\t" << gap << '\n';
	ofs << "End score:\t\t" << NWtable[xAxis - 1][yAxis - 1] << '\n';

	// print the calculated values from scorechecker and percentages
	ofs << "\nLength:\t\t\t" << vecAlignment[index0].size() << '\n';
	ofs << "Identity:\t\t" << iIdentityCount << "/"
		<< vecAlignment[index0].size() << " (" << identitypercent << "%)\n";
	ofs << "Similarity:\t\t" << iSimilarityCount + iIdentityCount << "/"
		<< vecAlignment[index0].size() << " (" << similaritypercent << "%)\n";
	ofs << "Gaps:\t\t\t" << iGapCount << "/"
		<< vecAlignment[index0].size() << " (" << gappercent << "%)\n";
	ofs << "#############################################################\n";
}



// writes alignment to file and flushes memory
void outputStream()
{
	std::cout << "Alignments made: " << iAlignment << '\r';
	// if the alignment is smaller than 60 pairs: do not use more lines
	if (vecAlignment[index0].size() < 60) {

		// first alignment sequence
		ofs << "\nAlignment " << iAlignment << '\n'
			<< vecAlignment[index0] << '\n';

		// check identity/similarity
		for (int i = 0; i < vecAlignment[index0].size(); ++i) {
			if (vecAlignment[index0][i] == vecAlignment[index1][i])
				ofs << "|";
			else if (bProtein == true)
				ofs << ConserveCheck
				(vecAlignment[index0][i], vecAlignment[index1][i]);
			else
				ofs << " ";
		}

		// output second alignment sequence
		ofs << "\n" << vecAlignment[index1] << std::endl;
	}

	// to make it legible for alignments over 60 pairs, use newlines
	else {
		ofs << "\nAlignment " << iAlignment << '\n';
		int i = 0; // counting variables: used later outside for-loops
		int j = 0;
		int k = 0;
		int subtract0 = 0; // to keep track of the sequence count
		int subtract1 = 0; // for both indices

		// print the first alignment
		for (; i < vecAlignment[index0].size(); ++i) {
			if (vecAlignment[index0][i] == '-')
				++subtract0;
			ofs << vecAlignment[index0][i];

			// when 60^n is reached, end the line 
			if ((i + 1) % 50 == 0 && i != 0) {
				ofs << " " << i - subtract0 + 1 << '\n';

				// output the match or mismatch pairs for 60^n
				for (k = j; k < i + 1; ++k) {
					if (vecAlignment[index0][k] == vecAlignment[index1][k])
						ofs << "|";
					else if (bProtein == true)
						ofs << ConserveCheck
						(vecAlignment[index0][k], vecAlignment[index1][k]);
					else
						ofs << " ";
				}
				ofs << "\n";

				// output the second sequence up to 60^n
				for (; j < i + 1; ++j) {
					if (vecAlignment[index1][j] == '-')
						++subtract1;
					ofs << vecAlignment[index1][j];
				}

				// new line and continue
				ofs << " " << i - subtract1 + 1 << "\n\n";
			}
		}
		ofs << " " << i - subtract0 << "\n";

		// print the remainder of the sequence and matches and flush
		for (; k < i; ++k) {
			if (vecAlignment[index0][k] == vecAlignment[index1][k])
				ofs << "|";
			else if (bProtein == true)
				ofs << ConserveCheck
				(vecAlignment[index0][k], vecAlignment[index1][k]);
			else
				ofs << " ";
		}
		ofs << '\n';

		for (; j < i; ++j) {
			if (vecAlignment[index1][j] == '-')
				++subtract1;
			ofs << vecAlignment[index1][j];
		}
		ofs << " " << j - subtract1 << std::endl;
	}

	// flush alignment from memory
	vecAlignment.erase(vecAlignment.begin(), vecAlignment.begin() + 2);
	++iAlignment; // alignment number
}



/************************************
*		 Alignment functions		*
*************************************/

// initialises the first alignment
void initialise()
{
	// start in the bottom right of the table
	x = vecSequences[seq0].size();
	y = vecSequences[seq1].size();
}



// appends the alignments based on gaps or matches
void alignmentMain(char& cDirection)
{
	if (cDirection == 'd') { // diagonal
		vecAlignment[index0].push_back(vecSequences[seq0][x-1]);
		vecAlignment[index1].push_back(vecSequences[seq1][y-1]);
		x -= 1; y -= 1;
	}
	else if (cDirection == 'l') { // left
		vecAlignment[index0].push_back(vecSequences[seq0][x-1]);
		vecAlignment[index1].push_back('-');
		x -= 1;
	}
	else { // up
		vecAlignment[index0].push_back('-');
		vecAlignment[index1].push_back(vecSequences[seq1][y-1]);
		y -= 1;
	}
}



// handles the alignment
void alignmentHandler()
{
	// iterate until end (top left) of matrix is reached, basically
	while (x > 0 || y > 0) { 

		// get the routes
		std::string routes;
		routes = alignmentTracker[x][y];

		// if there are multiple routes, there are alternative alignments
		while (routes.size() > 1) {
			if (altCount < maxAltCount) {
				++altCount;

				// store what has been aligned so far
				vecAlignment.push_back(vecAlignment[index0]);
				vecAlignment.push_back(vecAlignment[index1]);

				// store coordinate and alternative route
				altTracker.push_back({ x, y, routes.back() });
			}
			// remove the route
			routes = routes.substr(0, routes.size() - 1);
		}
		// only one route remains
		cDirection = routes[0];
		
		// align next
		alignmentMain(cDirection);
	}

	// check similarity, identity and gaps of the first alignment
	if (iAlignment == 1) {
		scoreChecker();
		outputHeader();
	}

	// alignments are done from back to front, so reverse them
	reverse(vecAlignment[index0].begin(), vecAlignment[index0].end());
	reverse(vecAlignment[index1].begin(), vecAlignment[index1].end());
	outputStream(); // output complete alignments to file

}



// handles all alternative alignments
void altHandler()
{
	while (!altTracker.empty()) {
		x = altTracker[0].x; // starting x
		y = altTracker[0].y; // starting y
		char cDirection = altTracker[0].route; // the route not yet taken
		alignmentMain(cDirection); // update x, y
		alignmentHandler(); // continue alignment normally
		altTracker.erase(altTracker.begin()); // delete the completed alt 
	} // continue until no more alts are left
}



/************************************
*			Main Function			*
*************************************/

int main(int argc, char* argv[])
{
	// get time at start
	const auto startTime = std::chrono::high_resolution_clock::now();

	// will not start if no data argument is given
	try {
		std::clog << "Starting program: " << argv[0] << '\n';
		if (argc == 1) {
			throw std::invalid_argument
			("No file name specified in arguments.\n");
		}

		// name and read the data
		std::clog << "\n\n********************\n" << "    Reading data\n"
			<< "********************\n";
		for (int i = 1; i < argc; ++i) {
			std::string dataFile = argv[i];
			dataReader(dataFile);
		}


		// cannot compare less than 2 sequences
		if (vecSequences.size() < 2)
			throw std::logic_error
			("less than 2 valid sequences found in the data file.\n");

		// output file stream
		std::ofstream ofs;
		ofs.open("AlignmentResult.txt");
		if (!ofs.is_open())
			throw std::runtime_error("unable to open output file.\n"
				"File may be in use by another programme.\n");

		// begins the alignment and outputs alignments
		initialise();
		std::clog << "\n\n********************\n" << "    Making table\n"
			<< "********************\n\n";
		tableMaker();
		std::clog << "\n********************\n" << " Aligning sequences\n"
			<< "********************\n\n";
		alignmentHandler();
		altHandler();

		// calculate time taken
		const auto endTime = std::chrono::high_resolution_clock::now();
		const auto duration =
			std::chrono::duration_cast<std::chrono::seconds>(endTime - startTime);

		// finalise
		if (yAxis < 10)
			showTable();
		std::clog << '\n' << "Alignment finished!\n";
		std::clog << "\nTotal time taken: " << duration.count() << " seconds.\n";
		ofs.close();

		if (ofs.is_open())
			throw std::runtime_error("unable to close output file.\n"
				"File may be in use by another programme.\n");
	}
	catch (std::exception& error) {
		std::cerr << "\nError: " << error.what();
		exit(EXIT_FAILURE);
	}

	return(0);
}