#!/usr/bin/perl
#
# Return the current 'tag' or branch name with no trailing newline.
#
# $Id$
#
use warnings;
use strict;
my $svnURL = '$HeadURL$';
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
