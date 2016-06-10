use Module::Build;
use strict;
use warnings;
use PDL::Core::Dev;

my $inc = PDL_INCLUDE;
$inc =~ s/^"?-I//;
$inc =~ s/"$//;

my $version = 0.02_01;

my $build = Module::Build->new(
	module_name => 'PDL::Fit::ExpRate',
	dist_version_from => 'lib/PDL/Fit/ExpRate.pm.PL',
	license  => 'perl',
	configure_requires => {
		# meta_merge was added in v 0.28:
		'Module::Build' => 0.28,
		PDL             => '2.4.0',
	},
	requires => {
		PDL  => '2.4.0',
		perl => '5.6.1',
	},
	needs_compiler => 1,
	include_dirs => $inc,
	meta_merge => {
		resources => {
			repository => 'http://github.com/run4flat/PDL-Fit-ExpRate',
			bugtracker => 'http://github.com/run4flat/PDL-Fit-ExpRate/issues',
		},
		no_index => {
			file => ['README.pod'],
			directory => [qw(examples t)],
		},
	},
	license => 'perl',
);

$build->create_build_script;
