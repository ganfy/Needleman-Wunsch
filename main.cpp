#include <iostream>
#include <string>
#include <queue>
#include <stack>
#include <vector>
#include <unordered_map>
#include <fstream>
#include <cctype>
#include <algorithm>
#include <climits>
using namespace std;

const int MAX_ALIGNMENTS = 1000; // Límite máximo de alineamientos a considerar

int countGaps(const string& seq) {
    int gapCount = 0;
    for (char c : seq) {
        if (c == '-') {
            gapCount++;
        }
    }
    return gapCount;
}

vector<pair<string, string>> Alignments;
string s1, s2, s3;
vector<vector<int>> score;
vector<vector<vector<short>>> paths;

/*Needleman and Wunsch*/
int NeedlemanWunsch() {
    int m = s1.length();
    int n = s2.length();

    score.assign(m + 1, vector<int>(n + 1, 0));
    paths.assign(m + 1, vector<vector<short>>(n + 1, vector<short>()));

    score[0][0] = 0;

    for (int i = 1; i <= m; i++){
      score[i][0] = score[i - 1][0] - 2;
      paths[i][0].push_back(1);
    }

    for (int j = 1; j <= n; j++){
       score[0][j] = score[0][j - 1] - 2;
       paths[0][j].push_back(2);
    }


    for (int i = 1; i <= m; i++) {
        for (int j = 1; j <= n; j++) {
            int score_a = score[i - 1][j - 1] + (s1[i - 1] == s2[j - 1]) * 2 - 1;
            int score_b = score[i - 1][j] - 2;
            int score_c = score[i][j - 1] - 2;
            int max_score = max(max(score_a, score_b), score_c);

            if (score_a == max_score)
                paths[i][j].push_back(0);
            if (score_b == max_score)
                paths[i][j].push_back(1);
            if (score_c == max_score)
                paths[i][j].push_back(2);

            score[i][j] = max_score;
        }
    }

    return score[m][n];
}

/*
void get_alignments() {
  int m = s1.length();
  int n = s2.length();

  stack<tuple<int, int, string, string>> stack;
  stack.push({m, n, "", ""});

  while (!stack.empty()) {
      auto i = get<0>(stack.top()), j = get<1>(stack.top());
    auto s1_align = get<2>(stack.top()), s2_align = get<3>(stack.top());
      stack.pop();

      if (i == 0 && j == 0) {
          reverse(s1_align.begin(), s1_align.end());
          reverse(s2_align.begin(), s2_align.end());

          Alignments.push_back({s1_align, s2_align});
      } else {

          for (auto direction : paths[i][j]) {
              string new_s1_align = s1_align;
              string new_s2_align = s2_align;

              if (direction == 0) {
                  new_s1_align += s1[i - 1];
                  new_s2_align += s2[j - 1];
                  stack.push({i - 1, j - 1, new_s1_align, new_s2_align});
              } else if (direction == 1) {
                  new_s1_align += s1[i - 1];
                  new_s2_align += '-';
                  stack.push({i - 1, j, new_s1_align, new_s2_align});
              } else if (direction == 2) {
                  new_s1_align += '-';
                  new_s2_align += s2[j - 1];
                  stack.push({i, j - 1, new_s1_align, new_s2_align});
              }
          }
      }
  }
}
*/

void get_alignments() {
    int m = s1.length();
    int n = s2.length();

    stack<tuple<int, int, string, string>> stack;
    stack.push({m, n, "",  ""});

    int minGapCount = INT_MAX;
    pair<string, string> bestAlignment; 

    int alignmentCount = 0, gapCount = INT_MAX;

    while (!stack.empty() && alignmentCount < MAX_ALIGNMENTS) {
      auto i = get<0>(stack.top()), j = get<1>(stack.top());
      auto s1_align = get<2>(stack.top()), s2_align = get<3>(stack.top());
        stack.pop();

        if (alignmentCount > 0 && i == 0 && j == 0)
          gapCount = countGaps(s1_align) + countGaps(s2_align);

        if (gapCount < minGapCount) {
            minGapCount = gapCount;
            bestAlignment = {s1_align, s2_align};
            //cout << minGapCount << " Best alignment: " << bestAlignment.first << endl;
        }

        if (++alignmentCount >= MAX_ALIGNMENTS) {
            break;
        }


      for (auto direction : paths[i][j]) {
          string new_s1_align = s1_align;
          string new_s2_align = s2_align;

          if (direction == 0) {
              new_s1_align += s1[i - 1];
              new_s2_align += s2[j - 1];
              stack.push({i - 1, j - 1, new_s1_align, new_s2_align});
          } else if (direction == 1) {
              new_s1_align += s1[i - 1];
              new_s2_align += '-';
              stack.push({i - 1, j, new_s1_align, new_s2_align});
          } else if (direction == 2) {
              new_s1_align += '-';
              new_s2_align += s2[j - 1];
              stack.push({i, j - 1, new_s1_align, new_s2_align});
          }
      }
    }


    cout << "Best Alignment (with minimum gaps):" << endl;
    cout << "Sequence 1: " << bestAlignment.first << endl;
    cout << "Sequence 2: " << bestAlignment.second << endl;
}


void get_alignments_R(string s1_align, string s2_align, int i, int j) {
  if (i == 0 && j == 0) {
      Alignments.push_back({s1_align, s2_align});
      return;
  }

  for (auto direction : paths[i][j]) {
      string new_s1_align = s1_align;
      string new_s2_align = s2_align;

      if (direction == 0) {
          new_s1_align = s1[i - 1] + new_s1_align;
          new_s2_align = s2[j - 1] + new_s2_align;
          get_alignments_R(new_s1_align, new_s2_align, i - 1, j - 1);
      } else if (direction == 1) {
          new_s1_align = s1[i - 1] + new_s1_align;
          new_s2_align = '-' + new_s2_align;
          get_alignments_R(new_s1_align, new_s2_align, i - 1, j);
      } else if (direction == 2) {
          new_s1_align = '-' + new_s1_align;
          new_s2_align = s2[j - 1] + new_s2_align;
          get_alignments_R(new_s1_align, new_s2_align, i, j - 1);
      }
  }
}

vector<string> read_sequences(const string& filename) {
    ifstream file(filename);
    if (!file.is_open()) {
        cerr << "Error: Unable to open file " << filename << endl;
        exit(1);
    }

    vector<string> sequences;
    string line, sequence;
    int virus = -1;
    while (getline(file, line)) {
        if (line.empty())
            continue;

        if (line.find_first_not_of(" \t\n\v\f\r") != string::npos && isalpha(line[0])) {
            if (!sequence.empty())
              sequences[virus] += sequence;
            sequence = ""; 
           sequences.push_back(sequence);
            virus++;
        } else {
          line.erase(remove_if(line.begin(), line.end(), ::isspace), line.end());
          line.erase(remove_if(line.begin(), line.end(), ::isdigit), line.end());
          sequence += line;
        }
    }
    if (!sequence.empty())
        sequences[virus] += sequence;

    file.close();
    return sequences;
}

void save_sequences(vector<string> seqs){
  ofstream file("sequences.txt");
  if (!file.is_open()){
    cerr << "Error: Unable to open file sequences.txt" << endl;
  }
  for (auto seq : seqs){
    file << seq << endl;
  }
}

int main() {
  s1 = "AAAC";
  s2 = "AGC";

  vector<string> sequences = read_sequences("Sequencias.txt"); 
  save_sequences(sequences);
  s1 = sequences[0].substr(0,1100);
  s2 = sequences[1].substr(0,1100);

  cout << s1 << endl << s2 << endl;

  int optimal_score = NeedlemanWunsch();
  cout << "Optimal alignment score: " << optimal_score << endl;

/*
for (auto dirs: paths){
  for (auto dir: dirs){
    for (auto d: dir)
      cout << d << ",";
    cout << "\t";
  }
  cout << endl;
}
*/

get_alignments();

  cout << "All possible optimal alignments:" << endl;
  cout << Alignments.size() << endl;

  /*
  for (const auto& alignment : Alignments) {
      cout << alignment.first << endl;
      cout << alignment.second << endl << endl;
  }
*/
}
