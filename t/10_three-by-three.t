# tests for three_by_three_Householder
use strict;
use warnings;
use Test::More tests => 11;
use PDL;
use PDL::Fit::ExpRate;

# Create a 3x3 diagonal matrix and make sure the we get the right answer
my $A = pdl q[1 0 0 ; 0 2 0 ; 0 0 3];
my $y = pdl q[5     ,   -4  ,     6];
my $coefs = three_by_three_Householder($A, $y);
my $expected = pdl(5, -2, 2);
ok(all ($coefs == $expected), 'Super-simple diagonalization works')
	or diag("Expected $expected but got $coefs");

# Build a random vector and matrix, multiply them, and then diagonalize
for (1..10) {
	$expected = zeroes(3)->grandom;
	$A = $A->grandom;
	$y = ($A x $expected->transpose)->squeeze;
	$coefs = three_by_three_Householder($A, $y);
	ok(all (approx($coefs, $expected)), 'Diagonalization works on random data')
		or diag("Expected $expected but got $coefs");
}
