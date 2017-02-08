#!/usr/bin/perl -w

######################################################################
########################################################################
## Copyright 2016 Fundação Oswaldo Cruz
## Author: Laura Rabelo Leite
## META_pipeline.pl is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by the Free Software Foundation, version 3 of the License.
##
## META_pipeline.pl is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY;
## without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
## See the GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License along with META_pipeline.pl (file: COPYING).
## If not, see <http://www.gnu.org/licenses/>.
##
########################################################################
########################################################################

use Tk;
use Tk::DirTree;
use Tk::Pane;
use Tk::LabFrame;
use Cwd;
use strict;

my $config;
my $input;
my $output;

# create the application window
my $mw = MainWindow->new(-screen=>$ARGV[0] || $ENV{DISPLAY} , -background => "GREY" );
$mw->Frame->pack(-side => 'top', -fill => 'x');


# add a window title
$mw->title("META_pipeline (2016)");

$mw->Button(-text => "Configuration file", -command => sub {
  my @types =
       (["Log files", [qw/.txt .log/]],
        ["All files",        '*'],
       );
  $config = $mw->getOpenFile(-filetypes => \@types) or return();
  warn "$config selected\n";
  })->
    pack(-side => 'right', -anchor => 'e');

$mw->Button(-text => "Input directory", -command => sub {
    $input = $mw->chooseDirectory(-initialdir => '~',
                                   -title => 'Choose a directory');
    if (!defined $input) {
        warn 'No directory selected';
    } else {
        warn "Selected $input";
    }
    })->
    pack(-side => 'right', -anchor => 'e');


$mw->Button(-text => "Output directory", -command => sub {
    $output = $mw->chooseDirectory(-initialdir => '~',
                                   -title => 'Choose a directory');
    if (!defined $output) {
        warn 'No directory selected';
    } else {
        warn "Selected $output";
    }
    })->
    pack(-side => 'right', -anchor => 'e');

$mw->Button(-text => "Exit", -command => sub { exit; } )->
    pack(-side => 'bottom', -anchor => 'e');

$mw->Button(-text => "Run", -command => sub { warn "Running: META_pipeline_CL.pl $input $config $output\n"; `META_pipeline_CL.pl $input $config $output`; } )->
    pack(-side => 'bottom', -anchor => 'e');

# wait until the user exits the app
MainLoop;

# exit the app
exit 0;

