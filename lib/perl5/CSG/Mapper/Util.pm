package CSG::Mapper::Util;

use base qw(Exporter);

use CSG::Base qw(file);
use CSG::Constants qw(:mapping);

our @EXPORT = ();

our @EXPORT_OK = (
  qw(
    parse_align_status
    parse_time
    )
);

our %EXPORT_TAGS = (
  parsers => [
    qw(
      parse_align_status
      parse_time
      )
  ],
);

sub parse_align_status {
  my $file = shift;

  my @results = ();
  for my $line (read_lines($file)) {
    chomp($line);
    my @parts = split(/\s/, $line);

    my $result_ref = {
      sampleid => $parts[0],
      orig_bam => $parts[1],
      results  => [],
    };

    if ($parts[2] =~ /\,/) {
      my @results = split(/,/, $parts[2]);
      my @states  = split(/,/, $parts[3]);

      for (0 .. $#results) {
        push @{$result_ref->{results}}, {result => $results[$_], state => $states[$_]};
      }
    } else {
      push @{$result_ref->{results}}, {result => $parts[2], state => $parts[3]};
    }

    push @results, $result_ref;
  }

  return @results;
}

sub parse_time {
  my ($time) = @_;
  return unless defined $time;

  for my $regexp (@TIME_FORMAT_REGEXPS) {
    return DateTime::Duration->new(%+) if $time =~ $regexp;
  }

  return undef;
}

1;
