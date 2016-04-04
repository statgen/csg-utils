package CSG::Constants;

use base qw(Exporter);
use CSG::Base;

our @EXPORT = (
  qw(
    $EMPTY
    $COMMA
    $UNDERSCORE
    $PERIOD
    $TRUE
    $FALSE
    $PIPE
    $DASH
    $SPACE
    $TIMEZONE
    )
);

our @EXPORT_OK = (
  qw(
    $EMPTY
    $COMMA
    $UNDERSCORE
    $PERIOD
    $TRUE
    $FALSE
    $PIPE
    $DASH
    $SPACE
    $MAX_DELAY
    $TIMEZONE
    $VALID_CLUSTER_REGEXPS
    @TIME_FORMAT_REGEXPS
    )
);

our %EXPORT_TAGS = (
  all => [
    qw(
      $EMPTY
      $COMMA
      $UNDERSCORE
      $PERIOD
      $TRUE
      $FALSE
      $PIPE
      $DASH
      $SPACE
      $TIMEZONE
      $VALID_CLUSTER_REGEXPS
      $MAX_DELAY
      @TIME_FORMAT_REGEXPS
      )
  ],
  basic => [
    qw(
      $EMPTY
      $COMMA
      $UNDERSCORE
      $PERIOD
      $TRUE
      $FALSE
      $PIPE
      $DASH
      $SPACE
      )
  ],
  mapping => [
    qw(
      $VALID_CLUSTER_REGEXPS
      $MAX_DELAY
      @TIME_FORMAT_REGEXPS
      )
  ],
);

Readonly::Scalar our $EMPTY      => q{};
Readonly::Scalar our $COMMA      => q{,};
Readonly::Scalar our $UNDERSCORE => q{_};
Readonly::Scalar our $PERIOD     => q{.};
Readonly::Scalar our $TRUE       => q{1};
Readonly::Scalar our $FALSE      => q{0};
Readonly::Scalar our $PIPE       => q{|};
Readonly::Scalar our $DASH       => q{-};
Readonly::Scalar our $SPACE      => q{ };
Readonly::Scalar our $TIMEZONE   => q{America/Detroit};

Readonly::Scalar our $MAX_DELAY             => 120;
Readonly::Scalar our $VALID_CLUSTER_REGEXPS => qr{csg|flux};

Readonly::Array our @TIME_FORMAT_REGEXPS => (

  # dd-hh:mm:ss or dd:hh:mm:ss
  qr/(?<days>\d{1,2})(?:\-|:)(?<hours>\d{2}):(?<minutes>\d{2}):(?<seconds>\d{2})/,

  # hhh:mm:ss
  qr/(?<hours>\d{1,3}):(?<minutes>\d{2}):(?<seconds>\d{2})/,

  # hh:mm
  qr/(?<hours>\d{1,2}):(?<minutes>\d{2})/,

  # sssssss
  qr/(?<seconds>\d{1,7})/,
);

1;
