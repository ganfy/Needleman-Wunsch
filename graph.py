import numpy as np
import matplotlib.pyplot as plt

def dot_matrix(string1, string2):
  matrix = np.zeros((len(string1), len(string2)))

  for i in range(len(string1)):
      for j in range(len(string2)):
          if string1[i] == string2[j]:
              matrix[i][j] = 1

  return matrix

def plot_dot_matrix(matrix, string1, string2):
  plt.imshow(matrix, cmap='coolwarm', interpolation='nearest')
  plt.xticks(np.arange(len(string2)), list(string2))
  plt.yticks(np.arange(len(string1)), list(string1))
  plt.xlabel('Influenza')
  plt.ylabel('SarsCov')
  plt.title('Dot Matrix of Strings')
  plt.colorbar(label='Matches')
  plt.show()


sequences = []

with open('sequences.txt', 'r') as file:
    for line in file:
        sequence = line.strip()
        sequences.append(sequence)
bacteria = sequences[0]
sarscov = sequences[1]
influenza = sequences[2]

#matrix = dot_matrix(bacteria, sarscov)
#plot_dot_matrix(matrix, bacteria, sarscov)
matrix = dot_matrix(sarscov, influenza)
plot_dot_matrix(matrix, sarscov, influenza)
#matrix = dot_matrix(sarscov, influenza)
#plot_dot_matrix(matrix, sarscov, influenza)
