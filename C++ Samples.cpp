
// Beispiel 1: Eine Funktion um ein Suffixarray zu konstruieren - um Laufzeit einzusparen, sortiert die Lambdafunktion ab Z.14 ff
// nur die numerischen Positionen, an denen das jeweilige Suffix im gegebenen Wort beginnt.
// Zur gleichen Aufgabe gibt es auch noch eine von mir programmierte MLR-Heuristik. 

void construct(std::vector<uint32_t>& sa, const std::string& text) { //sa = suffixarray
	sa.clear(); 
	
    //unsortiertes Befüllen des arrays
	if(!text.empty()){
		for (unsigned i = 0; i < text.length(); ++i) {
			sa.push_back(i); 
		}
    //Sortierung der Einträge
		sort(sa.begin(), sa.end(), [&text](int i, int j) {

			while (unsigned(i) < text.length() && unsigned(j) < text.length()) {

				if (text[i] < text[j])
				{
					return true;
				}
				else if (text[i] > text[j]) {

					return false;
				}
				else {
					++i;
					++j;
				}
			}
			if (i > j) {
				return true;
			}
			else {
				return false;
			}

		});
	}

//Beispiel 2: Matrizenbefüllung im Smith-Waterman-Algorithmus (für den Fall dass zwischen den zwei DNA Sequenzen, die verglichen werden, ein Match besteht)

//Initialisierung erste Spalte und erste Zeile: 
		for (unsigned i = 1; i <= s1.length(); i++) {
			fmatrix[0][i] = 0;
			tmatrix[0][i] = Traceback::HORIZONTAL;
		}	
	
		for (unsigned i = 1; i <= s2.length(); i++) {
			fmatrix[i][0] = 0;
			tmatrix[i][0] = Traceback::VERTICAL;
		}

	  //Initialisierung der restlichen Matrix:	
		for (unsigned i = 1; i <= s2.length(); i++){ //Abgehen erster Spalte 

			for (unsigned j = 1; j <= s1.length(); j++){ //Abgehen erster Zeile
			  // Match
				if (s1[j - 1] == s2[i - 1]) {  
				  //Fall 1: Diagonale + Matchscore ist höchste 
					if (fmatrix[i - 1][j - 1] + match >= fmatrix[i - 1][j] + gap
					 && fmatrix[i - 1][j - 1] + match >= fmatrix[i][j - 1] + gap
					 && fmatrix[i - 1][j - 1] + match >= 0) { 
						
						fmatrix[i][j] = fmatrix[i - 1][j - 1] + match;
						tmatrix[i][j] = Traceback::DIAGONAL;
					}

				  //Fall 3: Gap in S1
					else if (fmatrix[i][j - 1] + gap >= fmatrix[i - 1][j - 1] + match 
					 && fmatrix[i][j - 1] + gap >= fmatrix[i - 1][j] + gap
					 && fmatrix[i][j - 1] + gap >=	0) {

							fmatrix[i][j] = fmatrix[i][j - 1] + gap;
							tmatrix[i][j] = Traceback::HORIZONTAL;
					}
				  //Fall 4: Gap in S2
					else if (fmatrix[i - 1][j] + gap >= fmatrix[i - 1][j - 1] + match 
					 && fmatrix[i - 1][j] + gap >= fmatrix[i][j - 1] + gap
					 && fmatrix[i - 1][j] + gap >= 0) {
						
						fmatrix[i][j] = fmatrix[i - 1][j] + gap;
						tmatrix[i][j] = Traceback::VERTICAL;
					}
				  //Fall 5: 0 ist der höchste Wert
					else if (0 >= fmatrix[i - 1][j - 1] + match 
					 && 0 >= fmatrix[i - 1][j] + gap
					 && 0 >= fmatrix[i][j - 1] + gap) {
						 		
							fmatrix[i][j] = 0;
							tmatrix[i][j] = Traceback::NONE;
					}
				}

// Beispiel 3: Implementierung des BLAST-Algorithmus

    //Variablen: 
    vector<string> all_str;

 //Aufteilen der Query in word_size grosse q-Grame: 
    for (unsigned i = 0; i < query.length() - (word_size - 1); i++) {
        //Es wird nur bis zu dem Buchstaben iteriert, bei dem der Substring des Buchstabens noch wordsize lang ist. 

        all_str.push_back(query.substr(i, word_size)); // Die Lösungen werden in den Vector gepusht

    }

    vector<NHResult> my_vec;
    string word;
    vector<string> AS = { "A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y" };
    vector<string> list = { "" };
    int wscore = 0;
    vector<pair<string, int>> sorted;
    vector<string> temp;
    int i = 0;

    while (i < word_size) { // für jedes w-mer der gegebenen DNA-Sequenz werden alle Kombinationen erstellt

        temp.clear();
        for (unsigned k = 0; k < list.size(); k++) {

            for (unsigned c = 0; c < AS.size(); c++) {

                temp.push_back(list[k] + AS[c]);

            }
        }

        list = temp;
        i++;

    }

    for (unsigned x = 0; x < all_str.size(); x++) { //Überprüfung, ob jeweilige Neighborhoods den gegebenen Threshold knacken

        word = all_str[x];
        sorted.clear();

        for (unsigned l = 0; l < list.size(); l++) {

            string bet = list[l];
            wscore = 0;

            for (unsigned q = 0; q < word.size(); q++) {

                wscore += matrix.score(bet[q], word[q]); //die score-Funktion wurde vorgegeben, diese liest den score aus der gegebenen Substitutionsmatrix aus
            }

            if (wscore >= score_threshold) {

                sorted.emplace_back(bet, wscore); // Ergebnisse sollten in struct gefasst werden
            }

        }
        NHResult madeit{ word, sorted };
        my_vec.push_back(madeit);
    }

    return my_vec;
}

//Beispiel 4: main-Funktion, um Ergebnisse der Ausführung des BLAST Algorithmus und Laufzeit mittels knapper Eingabe direkt anzeigen zu lassen

int main(int argc, const char* argv[]) {

	if (argc != 6) {
		throw 1;
	}

	else {

		//für die Matrix
		ScoreMatrix m;
		m.load(argv[2]);

		//für die Zeitmessung
		double start;
		double end;

		int size = strtol(argv[3], NULL, 10);
		int threshold = strtol(argv[4], NULL, 10);
		int threads = strtol(argv[5], NULL, 10);

		if (threads != 1) {
			throw 1;
		}

		vector<NHResult> res;
		BLAST_Neighborhood bn;
        //Start der Zeitmessung
		start = omp_get_wtime();
		res = bn.generateNeighborhood(argv[1], m, size, threshold, threads);
        //Ende der Zeitmessung
		end = omp_get_wtime();

		// Ausgabe der Infixe mit Nachbarn und Scores in vorgegebenem Format
		for (unsigned i = 0; i < res.size(); i++) {

			cout << res[i].infix << ": ";

			for (unsigned j = 0; j < res[i].neighbors.size(); j++) {

				cout << "(" << res[i].neighbors[j].first << ", " << res[i].neighbors[j].second << ") ";
			}
			cout << endl;
		}

		// Ausgabe Zeit
		cout << "time: " << end - start << "s" << endl;

	}

}