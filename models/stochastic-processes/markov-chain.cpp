#include <iostream>
#include <vector>
#include <random>

using namespace std;

// Function to simulate a Markov Chain
vector<int> simulateMarkovChain(const vector<vector<double>>& transitionMatrix, int initialState, int steps, int seed = 42) {
    int numStates = transitionMatrix.size();
    vector<int> stateSequence;
    stateSequence.push_back(initialState);

    // Random number generator
    random_device rd;
    mt19937 generator(seed);
    uniform_real_distribution<double> distribution(0.0, 1.0);

    for (int i = 0; i < steps; ++i) {
        int currentState = stateSequence.back();
        double randNum = distribution(generator);
        double cumulativeProbability = 0.0;

        // Determine the next state based on transition probabilities
        for (int nextState = 0; nextState < numStates; ++nextState) {
            cumulativeProbability += transitionMatrix[currentState][nextState];
            if (randNum < cumulativeProbability) {
                stateSequence.push_back(nextState);
                break;
            }
        }
    }
    return stateSequence;
}

int main() {
    // Define states (Example: Weather States: 0 = Sunny, 1 = Cloudy, 2 = Rainy)
    vector<vector<double>> transitionMatrix = {
        {0.7, 0.2, 0.1}, // From Sunny
        {0.3, 0.4, 0.3}, // From Cloudy
        {0.2, 0.3, 0.5}  // From Rainy
    };

    int initialState = 0; // Start in 'Sunny' state
    int steps = 10;       // Number of transitions

    vector<int> stateSequence = simulateMarkovChain(transitionMatrix, initialState, steps);

    // Print the Markov Chain sequence
    cout << "Markov Chain Sequence (Weather States: 0=Sunny, 1=Cloudy, 2=Rainy):\n";
    for (int i = 0; i < stateSequence.size(); ++i) {
        cout << "Step " << i << ": State " << stateSequence[i] << "\n";
    }

    return 0;
}
