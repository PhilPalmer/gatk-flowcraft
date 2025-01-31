class Help {

    static def start_info(Map info, String time, String profile) {

        println ""
        println "============================================================"
        println "                F L O W C R A F T"
        println "============================================================"
        println "Built using flowcraft v1.4.2"
        println ""
        if (info.containsKey("fastq")){
        int nsamples = info.fastq / 2
        println " Input FastQ                 : $info.fastq"
        println " Input samples               : $nsamples"
        }
        if (info.containsKey("fasta")){
        println " Input Fasta                 : $info.fasta"
        }
        if (info.containsKey("accessions")){
        println " Input accessions            : $info.accessions"
        }
        println " Reports are found in        : ./reports"
        println " Results are found in        : ./results"
        println " Profile                     : $profile"
        println ""
        println "Starting pipeline at $time"
        println ""

    }

    static void complete_info(nextflow.script.WorkflowMetadata wf) {

        println ""
        println "Pipeline execution summary"
        println "=========================="
        println "Completed at                 : $wf.complete"
        println "Duration                     : $wf.duration"
        println "Success                      : $wf.success"
        println "Work directory               : $wf.workDir"
        println "Exit status                  : $wf.exitStatus"
        println ""

    }

    static def print_help(Map params) {

        println ""
        println "============================================================"
        println "                F L O W C R A F T"
        println "============================================================"
        println "Built using flowcraft v1.4.2"
        println ""
        println ""
        println "Usage: "
        println "    nextflow run main.nf"
        println ""
        println "       --fastq                     Path expression to paired-end fastq files. (default: $params.fastq) "
        println "       --reference                 Specifies the reference genome to be provided to HISAT2. (hisat2;haplotypecaller)"
        println "       --hisat2_index              Specifies the reference indexes to be provided to HISAT2. (hisat2)"
        println "       --hisat2_index_name         Specifies the reference indexes folder & basename to be provided to HISAT2, eg hisat2_index_folder/basename. (hisat2)"
        println "       --intervals                 Interval list file to specify the regions to call variants. (haplotypecaller)"
        
    }

}

class CollectInitialMetadata {

    public static void print_metadata(nextflow.script.WorkflowMetadata workflow){

        def treeDag = new File("${workflow.projectDir}/.treeDag.json").text
        def forkTree = new File("${workflow.projectDir}/.forkTree.json").text

        def metadataJson = "{'nfMetadata':{'scriptId':'${workflow.scriptId}',\
'scriptName':'${workflow.scriptName}',\
'profile':'${workflow.profile}',\
'container':'${workflow.container}',\
'containerEngine':'${workflow.containerEngine}',\
'commandLine':'${workflow.commandLine}',\
'runName':'${workflow.runName}',\
'sessionId':'${workflow.sessionId}',\
'projectDir':'${workflow.projectDir}',\
'launchDir':'${workflow.launchDir}',\
'startTime':'${workflow.start}',\
'dag':${treeDag},\
'forks':${forkTree}}}"

        def json = metadataJson.replaceAll("'", '"')

        def jsonFile = new File(".metadata.json")
        jsonFile.write json
    }
}