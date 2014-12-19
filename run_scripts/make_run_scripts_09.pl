use warnings;
use strict;

my $prescript = "09_perms_est_cross-tissue_h2";
open(QSUB, ">/group/im-lab/hwheeler/cross-tissue/qsub.txt");

for(my $i = 0; $i<=19; $i++){
    my $newi;
    if($i <= 9){
	$newi = 0 . $i;
    }else{
	$newi = $i;
    }
    my $outfile = "run_" . $prescript . "_subset" . $newi . ".sh";
    print QSUB "qsub run_scripts/$outfile\nsleep 2\n";
    open(OUT, ">$outfile");
    print OUT "#!/bin/bash\n#PBS -N R.gcta\n#PBS -S /bin/bash\n#PBS -l walltime=72:00:00\n#PBS -l nodes=1:ppn=4\n#PBS -l mem=4gb\n#PBS -o logs/\$\{PBS_JOBNAME\}.o\$\{PBS_JOBID\}.log\n#PBS -e logs/\$\{PBS_JOBNAME\}.e\$\{PBS_JOBID\}.err\n\ncd \$PBS_O_WORKDIR\nmodule load R\nmodule load gcta/1.24.4\n\ntime R --no-save < $prescript.r --args $newi\n";
}

