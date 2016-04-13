vector<int> KMP(string S) {
	int N = S.size();
	vector<int> F(N, 0);
	for (int i = 1; i < N; ++i) {
		int j = F[i];
		while (j && S[j] != S[i]) j = F[j];
		if (S[j] == S[i]) ++j;
		F[i + 1] = j;
    }
    return F;
}

void match(string S, vector<int> F, string X) {
	for (int i = 0, j = 0; i < (int) X.size(); ++i) {
		while (j && S[j] != X[i]) j = F[j];
		if (S[j] == X[i]) ++j;
		if (j == N) {
			// match at i - N + 1
			j = F[j];
		}
	}
}
