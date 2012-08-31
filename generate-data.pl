=head1 NAME

generate-data.pl - generate some exponential time series with additive noise

=head1 RUNNING

This script expects the following arguments:

 Switch    Meaning
 --A       additive constant
 --B       coefficient
 --tau     decay rate
 --noise   amplitude of the Gaussian noise
 --N       number of points to produce

A typical invocation of this script would look like this:

 generate-data.pl --A=150 --B=1 --tau=-5 --noise=2 --N=30

A file will be generated with a name constructed from the parameters. If this is
the seventh time you've run this, the output file will be called
F<A=150,B=1,N=30,noise=2,tau=-5,set=7.dat>.

=cut

use strict;
use warnings;

use Getopt::Long;
my	  ($A,  $B, $tau,   $amplitude, $N)
	= (150, 10, -10,    1,          30);

GetOptions(	"A=f"		=> \$A,
			"B=f"		=> \$B,
			"tau=f"		=> \$tau,
			"noise=f"	=> \$amplitude,
			"N=f"		=> \$N
		);

my $set_name = "A=$A,B=$B,N=$N,noise=$amplitude,tau=$tau";
my @set = glob("$set_name,set=*.dat");
my $out_file_name = "$set_name,set=" . (@set + 1) . '.dat';
#print "Creating file $out_file_name\n";

# Generate a bunch of data sets
use PDL;
my $xs = zeroes($N, 100)->xvals;
my $ys = $A + $B * exp($xs / $tau);# + $xs->grandom * $amplitude;

# Add some Gaussian noise
use PDL::GSL::RNG;
my $rng = PDL::GSL::RNG->new('taus');
$rng->set_seed(time());
$rng->ran_additive_gaussian($amplitude, $ys);

use PDL::Graphics::Prima::Simple;
#line_plot($xs, $ys);

#__END__

# Fit them to exponential trends
use PDL::Fit::ExpRate;
my @N_rounds_to_converge;
my @N_guarded_rounds;
my $guarded_rounds_counter = 0;
my ($As, $Bs, $taus) = fit_exp_rate($xs, $ys
	, threshold => 0.00001
	, iterations => 200
	#, trust_radius => 0.001
	, run_each_iteration => sub {
		my %args = @_;
		$guarded_rounds_counter++ if $args{forced_round};
		return 1;
	}
	, run_each_fit => sub{
		my %args = @_;
		push @N_rounds_to_converge, $args{N_rounds};
#		print "Fit $args{fit_count} (of $args{N_fits}) converged after $args{N_rounds} rounds\n";
#		print "$so_far / $total\n";
#		print "There were $guarded_rounds_counter guarded rounds\n";
		push @N_guarded_rounds, $guarded_rounds_counter;
		$guarded_rounds_counter = 0;
		return 1;
	}
);

use PDL::Graphics::Prima::Simple;

=pod

plot(
	-data => ds::Pair($xs, $ys),
	-fit  => ds::Func(sub {
			return $As + $Bs * exp( $_[0] / $taus );
		}
	),
	-real => ds::Func(sub {
			return $A + $B * exp( $_[0] / $tau );
		},
		colors => pdl(255, 0, 0)->rgb_to_color,
	),
	title => 'Fit',
);

__END__

=cut

my $N_rounds_to_converge = pdl(@N_rounds_to_converge);
plot(
	-data => ds::Pair($N_rounds_to_converge->hist($N_rounds_to_converge->minmax, 1),
		plotType => ppair::Histogram(),
	),
	x => { label => "Number of rounds to converge" },
	y => { label => 'Counts' },
);
my $N_guarded_rounds = pdl(@N_guarded_rounds);
plot(
	-data => ds::Pair($N_guarded_rounds->hist($N_guarded_rounds->minmax, 1),
		plotType => ppair::Histogram(),
	),
	x => { label => "Number of guarded rounds" },
	y => { label => 'Counts' },
);

plot(
	-data => ds::Pair($As->hist,
		plotType => ppair::Histogram(),
	),
	x => { label => "A value (should be $A)" },
	y => { label => 'Counts' },
);
plot(
	-data => ds::Pair($Bs->hist,
		plotType => ppair::Histogram(),
	),
	x => { label => "B value (should be $B)" },
	y => { label => 'Counts' },
);
plot(
	-data => ds::Pair($taus->hist,
		plotType => ppair::Histogram(),
	),
	x => { label => "tau value (should be $tau)" },
	y => { label => 'Counts' },
);

# Now let's plot the worst fit
my $worst = abs($Bs - $B)->maximum_ind;
my $worst_A = $As->at($worst);
my $worst_B = $Bs->at($worst);
my $worst_tau = $taus->at($worst);
my $worst_N = $N_rounds_to_converge->at($worst);
my $worst_guarded = $N_guarded_rounds->at($worst);
printf "Worst A: $worst_A; B: $worst_B; tau: $worst_tau; rounds: $worst_N; Guarded rounds: $worst_guarded\n";
use PDL::NiceSlice;
plot(
	-data => ds::Pair($xs(:, $worst;-), $ys(:, $worst;-)),
	-fit  => ds::Func(sub {
			return $worst_A + $worst_B * exp( $_[0] / $worst_tau );
		}
	),
	-real => ds::Func(sub {
			return $A + $B * exp( $_[0] / $tau );
		},
		colors => pdl(255, 0, 0)->rgb_to_color,
	),
	title => 'worst fit',
);

# idea: compute the DW statistic for the actual noise, and compare that with the
# DW statistic for the fit's residuals.
