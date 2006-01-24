package Geo::Coordinates::RDNAP;

use strict;

use vars qw($VERSION);
$VERSION = '0.01';

use Carp;
use Params::Validate qw/validate BOOLEAN SCALAR/;

use Exporter;
use vars qw/@ISA @EXPORT_OK/;
@ISA = qw/Exporter/;
@EXPORT_OK = qw/from_rd/;

my %a = (
    '01' => 3236.0331637,
    20 => -32.5915821,
    '02' => -0.2472814,
    21 => -0.8501341,
    '03' => -0.0655238,
    22 => -0.0171137,
    40 =>  0.0052771,
    23 => -0.0003859,
    41 =>  0.0003314,
    '04' =>  0.0000371,
    42 =>  0.0000143,
    24 => -0.0000090
);

my %b = (
    10 => 5261.3028966,
    11 => 105.9780241,
    12 =>  2.4576469,
    30 => -0.8192156,
    31 => -0.0560092,
    13 =>  0.0560089,
    32 => -0.0025614,
    14 =>  0.0012770,
    50 =>  0.0002574,
    33 => -0.0000973,
    51 =>  0.0000293,
    15 =>  0.0000291,
);

my %bessel = (
    a   => 6377397.155,
    e2  => 6674372e-9,
    f_i => 299.1528128
);

my %etrs89 = (
    a   => 6378137,
    e2  => 6694380e-9,
    f_i => 298.257222101
);

my %b2e = (
    tx  => 593.032,
    ty  => 26,
    tz  => 478.741,
    a   => 1.9848e-6,
    b   => -1.7439e-6,
    c   => 9.0587e-6,
    d   => 4.0772e-6
);

my @amersfoort_b = ( 3903_453.148, 368_135.313, 5012_970.306 );

sub from_rd {
    croak "Geo::Coordinates::RDNAP::from_rd needs two or three arguments"
        if (@_ !=2 && @_ != 3);

    my ($x, $y, $h) = (@_, 0);

    croak "Geo::Coordinates::RDNAP::from_rd: X out of bounds: $x"
        if ($x < -7 or $x > 300);
    croak "Geo::Coordinates::RDNAP::from_rd: Y out of bounds: $y"
        if ($y < 289 or $y > 629);

    # Use the approximated transformation.
    # Step 1: RD -> Bessel (spherical coords)

    $x = ($x/100) - 1.55;
    $y = ($y/100) - 4.63;

    my $lat = (52*60*60) + (9*60) + 22.178;
    my $lon = (5 *60*60) + (23*60) + 15.5;

    foreach my $i (keys %a) {
        my ($m, $n) = split //, $i;
        $lat += $a{$i} * ($x**$m) * ($y**$n);
    }

    foreach my $i (keys %b) {
        my ($m, $n) = split //, $i;
        $lon += $b{$i} * ($x**$m) * ($y**$n);
    }

    # Step 2: spherical coords -> X, Y, Z
    my @coords = _ellipsoid_to_cartesian($lat/3600, $lon/3600, $h, \%bessel);

    # Step 3: Bessel -> ETRS89
    @coords = _transform_datum( @coords, \%b2e, \@amersfoort_b );

    # Step 4: X, Y, Z -> spherical coords
    return _cartesian_to_ellipsoid(@coords, \%etrs89);
}

sub _to_rads {
    return $_[0] * 2*3.14159_26535_89793 /360;
}

sub _from_rads {
    return $_[0] / (2*3.14159_26535_89793) *360;
}

sub _ellipsoid_to_cartesian {
    my ($lat, $lon, $h, $para) = @_;

    my $sinphi = sin(_to_rads($lat));
    my $cosphi = cos(_to_rads($lat));
    my $n = $para->{a}/sqrt(1 - $para->{e2}*$sinphi*$sinphi);

    return (($n+$h)*$cosphi*cos(_to_rads($lon)),
            ($n+$h)*$cosphi*sin(_to_rads($lon)),
            ($n*(1-$para->{e2})+$h)*$sinphi );
}

sub _cartesian_to_ellipsoid {
    my ($x, $y, $z, $para) = @_;

    my $lon = atan2($y, $x);

    my $r = sqrt($x*$x+$y*$y);
    my $phi = 0;
    my $n_sinphi = $z;
    my $n;
    my $oldphi;

    do {
        $oldphi = $phi;
        $phi = atan2($z + $para->{e2}*$n_sinphi, $r);
        my $sinphi = sin($phi);
        $n = $para->{a}/sqrt(1-$para->{e2}*$sinphi*$sinphi);
        $n_sinphi = $n*$sinphi;
    } while (abs($oldphi-$phi) > 1e-8);

    my $h = $r/cos($phi) - $n;

    return (_from_rads($phi), _from_rads($lon), $h);
}

sub _transform_datum {
    my ($x, $y, $z, $t, $centre) = @_;

    return (
        $x + $t->{d}*($x-$centre->[0]) + $t->{c}*($y-$centre->[1])
            - $t->{b}*($z-$centre->[2]) + $t->{tx},
        $y - $t->{c}*($x-$centre->[0]) + $t->{d}*($y-$centre->[1])
            + $t->{a}*($z-$centre->[2]) + $t->{ty},
        $z + $t->{b}*($x-$centre->[0]) - $t->{a}*($y-$centre->[1])
            + $t->{d}*($z-$centre->[2]) + $t->{tz}
    );
}

1;
__END__

=head1 NAME

Geo::Coordinates::RDNAP - convert to/from Dutch RDNAP coordinate system

=head1 SYNOPSIS

  use Geo::Coordinates::RDNAP qw/from_rd/;

  # RD coordinates in km; height in meters
  my ($lat, $lon, $h) = from_rd( 150, 480, -2.75 );

  # lat/lon coordinates in degrees; height in meters
  # my ($x, $y, $h) = Geo::Coordinates::RD-NAP::to_rd( 52.75, 6.80, 10 );
  # (NOT YET IMPLEMENTED!)

=head1 DESCRIPTION

This module converts between two coordinate systems: RD-NAP and ETRS89.
ETRS89 is a geodesic frame of reference used in Europe, which is
approximately equal to the international reference frame WGS84.
GPS data. Coordinates in ETRS89 are given in degrees (latitude and
longitude) and meters (height above the reference ellipsoid).

RD-NAP is a Dutch coordinate system, consisting of the X and Y
coordinates of the Rijksdriehoekmeting, used e.g. in topographical maps,
and a Z coordinate which is the height above Normaal Amsterdams Peil,
the mean sea level at Amsterdam. X and Y are expressed in kilometers,
and Z in meters.

These transformations should only be used for locations in or close to
the Netherlands.

See http://www.rdnap.nl/ for a description of the RD-NAP system;
especially http://www.rdnap.nl/download/rdnaptrans.pdf for the formulas
used in this module.

=head2 Precision

This module implements an approximated transformation, which should be
accurate to about 25 cm in X and Y, and about 1 meter in the vertical
direction, for all locations in the Netherlands. The full
transformation, called RDNAPTRANS, is NOT implemented in this module. It
takes into account small deviations, measured at more than 5000 points
in the Netherlands.

Coordinates in ETRS89 deviate from WGS84 and ITRS because the former is
coupled to the Eurasian plate, which drifts a few cm per year compared
to other plates. The current (2006) difference between these coordinate
systems is in the order of 40 cm.

=head2 Disclaimer

Although this module implements conversion to/from the RD-NAP coordinate
system, it is not a product of RDNAP, the cooperation between the
Kadaster and Rijkwaterstaat, which maintains this coordinate system.

RDNAPTRANS is a trademark, presumably by Kadaster and/or
Rijkswaterstaat. This module is not an implementation of RDNAPTRANS. For
the official transformation software, visit http://www.rdnap.nl.

=head1 FUNCTIONS

=over 4

=item from_rd( $x, $y, $h )

Converts coordinates in the RD-NAP coordinate system to geographical
coordinates. The input are the X and Y coordinates in the RD system,
given in kilometers, and optionally the height above NAP in meters.

This should only be used for points in or close to the Netherlands. For
this area, X should roughly be between 0 and 300, and Y between 300 and 650.

The output is a list of three numbers: latitude and longitude in
degrees, according to the ETRS89 coordinate system, and height above the
ETRS89 reference geoid, in meters.

=back

=head1 BUGS

None known.

=head1 AUTHOR

Eugene van der Pijll C<< <pijll@cpan.org> >>

=head1 COPYRIGHT

Copyright (c) 2006 Eugene van der Pijll.
This program is free software; you can redistribute it and/or modify it
under the same terms as Perl itself.

The full text of the license can be found in the LICENSE file included
with this module.

=cut
