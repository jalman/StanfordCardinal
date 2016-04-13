public static BigInteger pollardRho(BigInteger y) {
	BigInteger x = new BigInteger(y.bitLength(), gen);
	BigInteger x2 = x;
	BigInteger c = new BigInteger(y.bitLength(), gen);
	BigInteger d;

	BigInteger two = BigInteger.valueOf(2);
	if (y.mod(two).equals(BigInteger.ZERO)) return two;

	int iters = 0;
	do {
		if (++iters > 3000) return BigInteger.ZERO;
		x = x.multiply(x).mod(y).add(c).mod(y);
		x2 = x2.multiply(x2).mod(y).add(c).mod(y);
		x2 = x2.multiply(x2).mod(y).add(c).mod(y);
		d = x.subtract(x2).gcd(y);
	} while (d.equals(BigInteger.ONE));

	return d;
}
