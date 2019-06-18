for f in assemblies/parts/*.fa; 
do 
	outdr=$(basename "$f" .fa)
	cat scripts/hpc_scripts/header.sb >scripts/hpc_scripts/annotation/$outdr.sb
    echo "prokka $f --force --outdir annotations/prokka/$outdr --norrna --notrna --metagenome --cpus 16" >>scripts/hpc_scripts/annotation/$outdr.sb
done
