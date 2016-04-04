package CSG::Storage::SlotCtl;

use App::Cmd::Setup -app;

sub global_opt_spec {
  return (
    ['debug|d',   'Debug output'],
    ['verbose|v', 'Verbose output'],
    ['dry-run',   'Dry run; show what would be done without actaully doing anything'],
  );
}

1;
