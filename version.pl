#!/usr/bin/perl
#
# Return the current 'tag' or branch name with no trailing newline.
#
# $Id$
#
use warnings;
use strict;
my $svnURL = '$HeadURL: svn+ssh://bw55@howdy.physics.nyu.edu/usr/local/svn/crange/branches/1.6/version.pl$';
if ($svnURL =~ m{/trunk/}) {
    my $svnVersion = qx{svnversion};
    chomp $svnVersion;
    $svnVersion =~ s{([0-9]+:)?([0-9]+)M?}{$2};
    print "trunk-r$svnVersion";
} else {
    my ($tag,$foo,$version) = ($svnURL =~ m{crange/(tag|branch)(s|es)/(.*)/});
    if ($version) {
        print $version;
    } else {
        print "NOTAG";
        exit 1;
    }
}
