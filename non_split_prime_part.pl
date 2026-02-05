sub aws {
    my ($V,$w,$ws) = @_;
    my $Ew = new Polytope(INEQUALITIES=>-$w|$V);
    my $Ews = new Polytope(INEQUALITIES=>-$ws|$V);
    my $P = new Polytope(POINTS=>ones_vector($V->rows)|$V);
    my $result = new Rational("inf");
    for my $x (@{rows($Ew->VERTICES)}) {
        my $result1 = new Rational("-inf");
        for my $x1 (@{rows($Ews->VERTICES)}) {
            my $result2 = new Rational("inf");
            my $idx2 = 0;
            for my $x2 (@{rows($P->VERTICES)}) {
                if ($x2 * $x1 != -$ws->[$idx2]) {
                    my $result3 = ($x2 * $x + $w->[$idx2]) / ($x2 * $x1 + $ws->[$idx2]);
                    if ($result3 < $result2) {
                        $result2 = $result3;
                    }
                }
                $idx2 = $idx2 + 1;
            }
            if ($result2 > $result1) {
                $result1 = $result2;
            }
        }
        if ($result1 < $result) {
            $result = $result1;
        }
    }
    return $result;
}

sub non_split_prime_part {
    my ($n, $w) = @_;
    my $h = hypersimplex(2,$n);
    my $s = splits($h->VERTICES, $h->GRAPH->ADJACENCY, $h->FACETS, $h->DIM);
    my $V = $h->VERTICES;
    my $sub = new fan::SubdivisionOfPoints(POINTS=>$V, WEIGHTS=>$w);
    my $i = 0;
    my $w_minus_w0 = zero_vector($w->dim);
    for my $j (@{splits_in_subdivision($V, $sub->MAXIMAL_CELLS, $s)}) {
        print $j;
        print "\n";
        my $row = $s->row($j);
        $ws = $V * $row;
        for my $i (0..$ws->dim-1) {
            if ($ws->[$i] < 0) {
                $ws->[$i] = 0;
            }
        }
        my $a = aws($V->minor(All, ~[0]), $w, $ws);
        $w_minus_w0 = $w_minus_w0 + $a * $ws;
        $i = $i + 1;
    }
    return $w_minus_w0;
}

$n = 5;
$w = new Vector([1,2,3,4,3,4,2,1,3/2,3/2]);
$w_minus_w0 = non_split_prime_part($n, $w);
print new fan::SubdivisionOfPoints(POINTS=>$V, WEIGHTS=>$w_minus_w0)->MAXIMAL_CELLS;