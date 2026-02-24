#!/usr/bin/env bash

set -euo pipefail

# Default/shared settings
num_children=100
num_replicates=3
include_indel=1
edit_distance=0
umic_threshold=0
thread=14
dispersion=2
tools_to_compare=UMI-nea,umi-tools,calib,humid,humid-fixed,humid-fixed-max

function performance() {
  # Set the variables
  pN=${num_founder}
  oN=${num_children}
  mut_ratio=${mutation_ratio}
  p1_ratio=${dispersion}
  do_indel=${include_indel}
  tools_to_run=${tools_to_compare}

  # Name of the current analysis
  name=sim_${pN}_${oN}_ul${umi_len}_err${err_rate}

  for rep in `seq 1 $num_replicates`; do
    truth=${name}/sim${rep}.truth.labels 
    if [ ! -f ${truth} ]; then
      echo "$truth MISSING" 1>&2
      continue
    fi

    # Simulated fout
    sim=${name}/sim${rep}.out 

    # Calculate common statistics for this run
    var_oN=`echo "$oN*$p1_ratio" | bc`
    total_umi=`cat ${sim}| wc -l`
    s_mu=`cat ${truth} | awk '{print $2}' | sort | uniq -c | awk '{print $1}' | awk -v n=0 '{n+=$1}END{print n/NR}'`
    s_var=`cat ${truth} | awk '{print $2}' | sort | uniq -c | awk '{print $1}' | awk -v a=0 -v m=$s_mu '{s=($1-m)*($1-m);a+=s}END{print int(a/(NR-1))}'`
    umi_err=`cat ${sim} | awk '$3!=""' | wc -l`
    umi_sub=`cat ${sim} | awk '$3!=""{gsub("[0-9]","",$3);gsub(":","",$3);split($3,n,"_");idl=0;for(i in n){m=n[i];split(m,n1,"/");if(length(n1[2])!=1){idl=1;break}};if(idl==0){print}}' | wc -l`
    umi_idl=`echo "$umi_err-$umi_sub" | bc`
    elist=`cat ${sim} | awk -v s=0 -v i=0 -v d=0 '$3!=""{gsub("[0-9]","",$3);gsub(":","",$3);split($3,n,"_");for(x in n){m=n[x];split(m,n1,"/");if(length(n1[2])==1){s+=1}else if(length(n1[2])>1){i+=1}else{d+=1}}}END{print i,d,s}'`
    bp_ins=`echo $elist | awk '{print $1}'`
    bp_del=`echo $elist | awk '{print $2}'`
    bp_sub=`echo $elist | awk '{print $3}'`
    bp_idl=`echo "$bp_ins+$bp_del" | bc`
    bp_ratio=`echo "$bp_ins $bp_del $bp_sub" | awk -v ind=$do_indel '{if(ind>0){a=$1;for(i=1;i<=3;i++){a=($i<a?$i:a)};printf "%3.3f-",$1/a;printf "%3.3f-",$2/a;printf "%3.3f",$3/a}else{print "NA"}}'`
    uniq_umi=`cat ${sim} | cut -f1 | sort | uniq | wc -l`

    for eval_t in `echo $tools_to_run | sed 's/,/\n/g'`; do
      # Set the number of threads
      if [ $eval_t == "umi-tools" ] || [ $eval_t == "humid" ] || [ $eval_t == "humid-fixed" ] || [ $eval_t == "humid-fixed-max" ]; then
          td=1
      else
          td=$thread
      fi

      # Output files for this tool/rep/etc
      scoref=${name}/$eval_t.sim${rep}.t$td.score
      tool_labels=${name}/$eval_t/sim${rep}.t$td.labels
      time=${name}/$eval_t.time
      log=${name}/log/$eval_t.sim${rep}.t$td.log

      if [ ! -f ${tool_labels} ];then
        echo "${tool_labels} MISSING" 1>&2
        continue
      fi

      # echo $scoref, ${tool_labels}, ${time}, ${log}
      # Check that this replicate is present in the time file
      
      set +e
      grep "$name r=$rep t=$td" ${time}
      exit_code=$?
      set -e

      if [ $exit_code -eq 1 ];then
        echo "rep $rep not found in $time" 1>&2
        continue
      fi

      # Calculate stats
      runtime_t="NA"
      rpu_cutoff="NA"
      rpu_model="NA"
      score_v=`if [ -f $scoref ]; then cat $scoref | awk '{printf "%.4f\n",$NF}'; else echo "NA"; fi`
      score_h=`if [ -f $scoref ]; then cat $scoref | awk '{printf "%.4f\n",$4}'; else echo "NA"; fi`
      score_c=`if [ -f $scoref ]; then cat $scoref | awk '{printf "%.4f\n",$7}'; else echo "NA"; fi`
      n_cluster=`cat ${tool_labels} | cut -d" " -f2 | sort -u | wc -l`
      est_mol=$n_cluster

      case $eval_t in
      UMIC-seq)
          runtime_t=`cat ${time} | grep -A 6 "$name r=$rep t=$td" | awk 'NR==3 || NR==7{print $2}' | awk -v n=0 '{split($NF,a,"m");n+=a[1]*60;split(a[2],b,"s");n+=b[1]}END{printf "%.2f\n",n}'`
          maxdist=`cat ${log} | head -1 | awk '{print $3}'`
          ;;
      calib)
          runtime_t=`cat ${time} | grep -A 2 "$name r=$rep t=$td" | tail -1 | awk '{split($NF,a,"m");n+=a[1]*60;split(a[2],b,"s");n+=b[1];printf "%.2f\n",n}'`
          maxdist=`cat ${log} | grep "error_tolerance" | head -1 | awk '{print $2}'`
          ;;
      umi-tools)
          runtime_t=`cat ${time} | grep -A 2 "$name r=$rep t=$td" | tail -1 | awk '{split($NF,a,"m");n+=a[1]*60;split(a[2],b,"s");n+=b[1];printf "%.2f\n",n}'`
          maxdist=`cat ${log} | grep "# threshold" | awk '{print $4}'`
          ;;
      UMI-nea)
          runtime_t=`cat ${time} | grep -A 2 "$name r=$rep t=$td" | tail -1 | awk '{split($NF,a,"m");n+=a[1]*60;split(a[2],b,"s");n+=b[1];printf "%.2f\n",n}'`
          maxdist=`cat ${log} | grep "maxdist" | awk '{print $NF}'`
          rpu_cutoff=`cat ${name}/$eval_t/sim${rep}.t$td.clustered.estimate | grep "rpu_cutoff" | awk '{print $NF}'`
          rpu_model=`cat ${name}/$eval_t/sim${rep}.t$td.clustered.estimate | head -1 | awk '{print ($1~/NB/?"negbinom":"kneeplot")}'`
          est_mol=`cat ${name}/$eval_t/sim${rep}.t$td.clustered.estimate | grep "estimated_molecules" | awk '{print $NF}'`
          ;;
      humid)
          runtime_t=`cat ${time}| grep -A 2 "$name r=$rep t=$td" | tail -1 | awk '{split($NF,a,"m");n+=a[1]*60;split(a[2],b,"s");n+=b[1];printf "%.2f\n",n}'`
          maxdist=`cat ${log} | grep "Edit distance" | awk '{print $3}'`
          ;;
      humid-fixed)
          runtime_t=`cat ${time}| grep -A 2 "$name r=$rep t=$td" | tail -1 | awk '{split($NF,a,"m");n+=a[1]*60;split(a[2],b,"s");n+=b[1];printf "%.2f\n",n}'`
          maxdist=`cat ${log} | grep "Edit distance" | awk '{print $3}'`
          ;;
      humid-fixed-max)
          runtime_t=`cat ${time}| grep -A 2 "$name r=$rep t=$td" | tail -1 | awk '{split($NF,a,"m");n+=a[1]*60;split(a[2],b,"s");n+=b[1];printf "%.2f\n",n}'`
          maxdist=`cat ${log} | grep "Edit distance" | awk '{print $3}'`
          ;;
      *)
          echo "invalid tool"
      esac
      echo "$pN $oN $var_oN $umi_len $err_rate $mut_ratio $rep $total_umi $s_mu $s_var $bp_sub $bp_idl $bp_ratio $umi_sub $umi_idl $uniq_umi $eval_t $maxdist $td $runtime_t $n_cluster $score_v $score_h $score_c $rpu_cutoff $rpu_model $est_mol"
    done
  done

  # echo "$pN $oN $var_oN $umi_len $err_rate $mut_ratio $rep $total_umi $s_mu $s_var $bp_sub $bp_idl $bp_ratio $umi_sub $umi_idl $uniq_umi $eval_t $maxdist $td $runtime_t $n_cluster $score_v $score_h $score_c $rpu_cutoff $rpu_model $est_mol"

}

# First, we print the header
  echo "num_founder mean_children_num variance_children_num umi_len err_rate insertion-deletion-substitution replicate total_umi simulated_mean_children_num simulated_variance_children_num substitution_base indel_base simulated_insertion-deletion-substitution substitution_only_umi indel_umi uniq_umi tool clustering_threshold thread runtime_in_sec dedup_umi_cluster V-measure homogeneity_score completeness_score RPU_cutoff RPU_cutoff_model estimated_molecule"

# TESTING
# num_founder=10000
# umi_len=50
# err_rate=0.01
# mutation_ratio=1-1-1
# tools_to_compare=humid-fixed-max
#
# performance
#
# exit
# END TESTING
#
# Run the Illumina benchmarks
mutation_ratio=1-1-40
for err_rate in 0.005 0.001; do
  for umi_len in 12 18; do
    for num_founder in 1000 10000; do
      performance
    done
  done
done

# Run the PacBio/ONT benchmarks
mutation_ratio=1-1-1
for err_rate in 0.01 0.03; do
  for umi_len in 25 50; do
    for num_founder in 1000 10000; do
      performance
    done
  done
done
