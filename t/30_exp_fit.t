use strict;
use warnings;
use Test::More;
use PDL;
use PDL::Fit::ExpRate;

##################
# Small data set #
##################
my $xs = sequence(30) + 100;
my $ys = 150 + 10 * exp($xs / -10);

my ($As, $Bs, $taus) = fit_exp_rate($xs, $ys
	, threshold => 0.00001
	, iterations => 200
	, trust_radius => 0.1
#	, run_each_iteration => sub {
#		my %info = @_;
#		return if $info{round} % 5 != 0;
#		print "Round $info{round}:\n";
#		delete $info{round};
#		for my $k (sort keys %info) {
#			print "  $k: $info{$k}\n";
#		}
#	},
);
my $expected = pdl(150, 10, -10);
my $coefs = cat($As, $Bs, $taus)->flat;
if (0 == $coefs->nbad) {
	ok(all (approx($coefs, $expected, 1e-2)), 'Exponential fitting on small noise-less dataset works')
		or do {
			diag("Got $coefs but expected $expected");
			($As, $Bs, $taus, my $err) = exp_fit_estimate_parameters($xs, $ys);
			diag("Quadratic estimate gave A=$As, B=$Bs, and lambda=$taus; err=$err");
		};
}
else {
	fail('Exponential fitting on small noise-less dataset works');
	diag "  -> Failed to converge";
}

###################################
# Large data set with minor noise #
###################################
$xs = sequence(30000);
my $tau = -1e5;
$ys = 150 + 10 * exp($xs / $tau) + $xs->grandom*0.001;
($As, $Bs, $taus) = fit_exp_rate($xs, $ys
	, threshold => 0.0001
	, iterations => 400
	, trust_radius => 0.1,
#	, run_each_iteration => sub {
#		my %info = @_;
#		return if $info{round} % 5 != 0;
#		print "Round $info{round}:\n";
#		delete $info{round};
#		for my $k (sort keys %info) {
#			print "  $k: $info{$k}\n";
#		}
#	},
);
$coefs = cat($As, $Bs, $taus)->flat;
if (0 == $coefs->nbad) {
	$expected = pdl(150, 10, $tau);
	ok(all (abs(($coefs - $expected) / $expected) < 1e-3)
		, 'Exponential fitting with small noise and lots of data works')
		or do {
			diag("Got $coefs but expected $expected");
			($As, $Bs, $taus, my $err) = exp_fit_estimate_parameters($xs, $ys);
			diag("Quadratic estimate gave A=$As, B=$Bs, and lambda=$taus; err=$err");
		};
}
else {
	fail('Exponential fitting with small noise and lots of data works');
	diag "  -> Failed to converge";
}

################################
# Large data set with no noise #
################################
$ys = 150 + 10 * exp($xs / $tau);
($As, $Bs, $taus) = fit_exp_rate($xs, $ys
	, threshold => 0.0001
	, iterations => 400
	, trust_radius => 0.1
#	, run_each_iteration => sub {
#		my %info = @_;
#		return if $info{round} % 5 != 0;
#		print "Round $info{round}:\n";
#		delete $info{round};
#		for my $k (sort keys %info) {
#			print "  $k: $info{$k}\n";
#		}
#	},
);
$coefs = cat($As, $Bs, $taus)->flat;
$expected = pdl(150, 10, $tau);
if (0 == $coefs->nbad) {
	ok(all (abs(($coefs - $expected) / $expected) < 1e-3)
		, 'Exponential fitting with no noise and lots of data works')
		or do {
			diag("Got $coefs but expected $expected");
			($As, $Bs, $taus, my $err) = exp_fit_estimate_parameters($xs, $ys);
			$taus = 1/$taus;
			diag("Quadratic estimate gave A=$As, B=$Bs, and tau=$taus; err=$err");
		};
}
else {
	fail('Exponential fitting with no noise and lots of data works');
	diag "  -> Failed to converge";
}

done_testing;
