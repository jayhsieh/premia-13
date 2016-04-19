#!/usr/bin/perl -w


#########################################################################
# Written and (C) by Jérôme Lelong <jerome.lelong@gmail.com>            #
#                                                                       #
# This program is free software; you can redistribute it and/or modify  #
# it under the terms of the GNU General Public License as published by  #
# the Free Software Foundation; either version 3 of the License, or     #
# (at your option) any later version.                                   #
#                                                                       #
# This program is distributed in the hope that it will be useful,       #
# but WITHOUT ANY WARRANTY; without even the implied warranty of        #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         #
# GNU General Public License for more details.                          #
#                                                                       #
# You should have received a copy of the GNU General Public License     #
# along with this program.  If not, see <http://www.gnu.org/licenses/>. #
#########################################################################

##
## This script finds the .tst files and runs the compute[] method on
## the Premia object found in the .tst file. The newly comuted valued
## is compared to the old stored in the .tst file. At the end of the
## tests, a mail is sent with the list of failures.
##
## This script must be run from the directory scripts !!
## Run perl regression_status.pl --help to print the help.

use strict;

use File::Find;
use File::Basename;
use Mail::Sendmail;
use Getopt::Long;


# for the convenience of &wanted calls, including -eval statements:
use vars qw/*name *dir *prune/;
*name   = *File::Find::name;
*dir    = *File::Find::dir;
*prune  = *File::Find::prune;

my $update=0;
my $test_dir='';
my $help=0;
my $email_addrs = 'jerome.lelong@imag.fr';
my $no_mail = 0;
my $ntests = 0;
my $nfail = 0;
my @failures = ('');
my $diff_file = "diff";
my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);

GetOptions(
           "update" => \$update,
           "help" => \$help,
           "email=s" => \$email_addrs,
           "no-mail" => \$no_mail,
           "test-dir=s" => \$test_dir,
          );

$test_dir =  '../tests/' unless ($test_dir);
$year = $year + 1900;
$mon = $mon + 1;
my $date = sprintf ("%4d%02d%02d%02d%02d", $year, $mon, $mday, $hour, $min);
$diff_file = join ('-', $diff_file, $date) . '.txt';


my %mail = (
            To => $email_addrs,
            From => 'viso@wizoo.inria.fr',
            Smtp => 'smtp.inria.fr', 
           );


sub dummyprint { print @_; }


##
## Prints the manual page
##
sub print_help {
    print("usage : regression_status.pl [--help] [--update] [--email=addr]
                                        [--test-dir=path]
\t--help : prints this help.
\t--email : specifies who receives the mail. addr must a
\t          comma separated list of email addresses.
\t--no-mail : do not send any email.
\t--update : performs a svn update and recompiles Premia
\t           before running the tests.
\t--test-dir : if specified, only checks the tests available 
\t             in path. path is relative.\n");
}

##
## Updates the svn, cleans everything and recompiles all
## The return value of every command is checked.
##
sub remake {
    my $makenfail = 0;
    my $makefailmsg = '';
    my $commandsvn="svn up";
    my $commandgen="./autogen.sh";
    my $commandclean="make distclean";
    my $commandmake="make install";
    chdir("../../");
    print "Running $commandsvn \n";
    if ($makenfail == 0 && system($commandsvn) != 0) {
        $makefailmsg = " error in svn";
        $makenfail ++;
    }
    print "Running $commandclean \n";
    if ($makenfail == 0 && system($commandclean)!= 0) {
        $makefailmsg = " error in clean";
        $makenfail ++;
    }
    print "Running $commandgen \n";
    if ($makenfail == 0 && system($commandgen)!= 0) {
        $makefailmsg = " error in autogen";
        $makenfail ++;
    }
    print "Running $commandmake \n";
    if ($makenfail == 0 && system($commandmake)!= 0) {
        $makefailmsg = " error in make";
        $makenfail ++;
    }
    chdir("nsp/scripts/");
    return (($makenfail, $makefailmsg)); 
}

sub wanted {
    if (/^.*\.tst\z/s) {
        $ntests ++;
        my @commands = ();
        my $command;
        push (@commands, "nsp -nw -e ");
        push (@commands,
              "\"exec(\\\"../libpremia/loader.sce\\\");  premia_init(); exec(\\\"compare_eps.sci\\\");"
             );
        push (@commands, ("fic = \\\"", $File::Find::name, "\\\"; load(fic);",
                          "Lold = P.get_method_results[];",
                          "P.compute[]; Lnew = P.get_method_results[];",
                          "if compare_eps(Lold,Lnew,\\\"", $diff_file, "\\\", fic)==0;",
                          "then printf(\\\"Status1\\\") end; quit;\"")
             );
        $command = join('', @commands);
        dummyprint (basename ($File::Find::name) . "...");
        my $code = ` $command `;
        if ($code =~ m/Status1/) {
            dummyprint (" OK \n");
        } else {
            dummyprint (" FAIL \n");
            push (@failures, basename ($File::Find::name));
            $nfail ++;
        }
    }
}

##
## Find and test all tst files. Compose a mail with the result of the
## tests.
##
sub test_all {
    my %find_opt;
    # this option is very important otherwise inside the function wanted,
    # we are chdired to the dirname of the current file.
    $find_opt{"no_chdir"} = 1;
    $find_opt{"wanted"} = \&wanted;

    File::Find::find(\%find_opt, $test_dir);
    my $message = join('',($nfail, " failure(s) out of ", $ntests, " tests\n"));
    $message = join( "\n", ($message, @failures));

    $mail{'Subject'} = "[Premia Tests]";
    $mail{'Message'} = $message;

    print ("\n******************\n");
    print ("** Test Summary **\n");
    print ("******************\n");
    print ($message);
}

sub do_sendmail {
    if (sendmail (%mail)) {
        print "Mail sent OK.\n";
    } else {
        print "Error sending mail: $Mail::Sendmail::error \n";
    }
}

if ($help == 1){
    print_help ();
    exit;
}

if ($update == 1) {
    my ($nfail, $msg) = remake();
    if ($nfail > 0) {
        $mail{'Subject'} = "[Premia Tests : ERROR IN COMPILATION !!!]";
        $mail{'Message'} = $msg;
        do_sendmail ();
        exit;
    }
}
test_all();
do_sendmail() unless ($no_mail == 1);

exit;

