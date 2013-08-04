#!/usr/bin/perl


# alarm-binary2: solves perfectly after 2K data points. .01
# alarm-binary: solves perfectly after 2K data points. .01
# alarm-binary-full: even with 10K still has some errors. It's because epsilon is too small. Need to recreate. One can sample CPDs such that epsilon is guaranteed to be non-trivial, by sampling marginals uniform distribution, and then sampling t from the appropriate range. However, epsilon will actually be smaller because when conditioning on parents/children it will give more information.
#
# The problem is that, because of the way we sampled the distributions, BIC does very badly. Likelihood does not suffice to distinguish the true network from the old network -- not a p-map. 

# compute PDAG first of the two graphs, then run Algorithm 4 of MMHC.pdf.

#$data = '../data/synthetic_examples/alarm-binary.dat';
#$outputdir = 'tmp_BIC';

$datadir = '../data/synthetic_examples/experiments/';
$outputdir = '../results/experiments/';
$method = 'sparsity_eps005';

for($exp=0; $exp<=9; $exp++) {

	print "\n\nEXPERIMENT #$exp\n\n";

	$datadir_full = $datadir . $exp . "/";
	$outputdir_partial = $outputdir . $exp;
	$outputdir_full = $outputdir . $exp . "/" . $method;
	if(!(-d $outputdir_partial)) {
		mkdir($outputdir_partial);
	}
	if(!(-d $outputdir_full)) {
		mkdir($outputdir_full);
	}

	$N_start = 4400;
	$N_end = 8000;
	$N_increment = 400;

	for($N=$N_start; $N<=$N_end; $N+=$N_increment) {

		print "\n$N\n";
		
		$time_file = $outputdir_full . "/bscore_time$N.dat";
		$model_file = $outputdir_full . "/model$N.mod";
		$mec_file = $outputdir_full . "/output$N.txt";
		$mec_file2 = $outputdir_full . "/bn$N.mec";

		$data_file = $datadir_full . "alarm" . $N . ".dat";

		# Create model file from this
		$start_time = time;
		#    `./c2/bscore -bic -nodes 37 -maxpa 4 -data $data_file -mod_out $model_file`;
		#    #    `./c2/bscore -alpha 10 -nodes 37 -maxpa 4 -data $data_file -mod_out $model_file`;
		`./c2/bscore -bic -nodes 37 -maxpa 4 -data $data_file -mod_out $model_file -edge_scores .005`;
		#    `./c2/bscore -bic -nodes 37 -maxpa 4 -data $data_file -mod_out $model_file -edge_scores .005`;
		print "./c2/bscore -bic -nodes 37 -maxpa 4 -data $data_file -mod_out $model_file -edge_scores .02\n";
		print " ";

		$end_time = time;
		$total_time = $end_time - $start_time;
		open(FILE, ">$time_file");
		print FILE $total_time;
		close(FILE);

		# Create settings file
		`cp gobnilp.set $outputdir_full/gobnilp$N.set`;
		`sed -i '' 's|%DIR|$outputdir_full|' $outputdir_full/gobnilp$N.set`;  #in-place editing: '' is the SUFFIX in -i[SUFFIX]
		`sed -i '' 's/%N/$N/' $outputdir_full/gobnilp$N.set`;

		# Run (note, need "gobnilp.set") in same directory
		print "gobnilp1.3/bin/gobnilp -g$outputdir_full/gobnilp$N.set $model_file > $mec_file";
		`gobnilp1.3/bin/gobnilp -g$outputdir_full/gobnilp$N.set $model_file > $mec_file`;

	}
}
