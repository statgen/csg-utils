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
    $TAB
    $COLON
    $SLASH
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
    $TAB
    $COLON
    $SLASH
    $MAX_DELAY
    $TIMEZONE
    $VALID_CLUSTER_REGEXPS
    $FASTQ_SUFFIX
    %CLUSTER_MAP
    @TIME_FORMAT_REGEXPS
    @READ_GROUP_FIELDS
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
      $TAB
      $COLON
      $SLASH
      $TIMEZONE
      $VALID_CLUSTER_REGEXPS
      $FASTQ_SUFFIX
      %CLUSTER_MAP
      $MAX_DELAY
      @TIME_FORMAT_REGEXPS
      @READ_GROUP_FIELDS
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
      $TAB
      $COLON
      $SLASH
      )
  ],
  mapping => [
    qw(
      $VALID_CLUSTER_REGEXPS
      $FASTQ_SUFFIX
      %CLUSTER_MAP
      $MAX_DELAY
      @TIME_FORMAT_REGEXPS
      @READ_GROUP_FIELDS
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
Readonly::Scalar our $TAB        => qq{\t};
Readonly::Scalar our $COLON      => q{:};
Readonly::Scalar our $SLASH      => q{/};
Readonly::Scalar our $TIMEZONE   => q{America/Detroit};

Readonly::Scalar our $MAX_DELAY             => 120;
Readonly::Scalar our $VALID_CLUSTER_REGEXPS => qr{csg|flux|dummy};
Readonly::Scalar our $FASTQ_SUFFIX          => q{.fastq.gz};

Readonly::Hash our %CLUSTER_MAP => (
  'red hat' => 'flux',
  'ubuntu'  => 'csg',
);

Readonly::Array our @READ_GROUP_FIELDS   => (qw(ID PL PU LB DS DT SM CN));
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
