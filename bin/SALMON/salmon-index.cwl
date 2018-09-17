cwlVersion: v1.0
class: CommandLineTool

requirements:
- class: InlineJavascriptRequirement

hints:
- class: DockerRequirement
  dockerPull: combinelab/salmon

inputs:

  fasta_file:
    type:
      - File
    inputBinding:
      position: 5
      prefix: "-t"

  kmer_len:
    type:
      - int?
    inputBinding:
      position: 6
      prefix: "-k"

  output_filename:
    type:
      - string?
    inputBinding:
      position: 7
      prefix: "-i"

  gencode:
    type:
      - boolean?
    inputBinding:
      position: 8
      prefix: "--gencode"

  keep_duplicates:
    type:
      - boolean?
    inputBinding:
      position: 9
      prefix: "--keepDuplicates"

  perfect_hash:
    type:
      - boolean?
    inputBinding:
      position: 10
      prefix: "--perfectHash"

  output_type:
    type:
      - string?
    inputBinding:
      position: 5
      prefix: "--type"

  sasamp:
    type:
      - int?
    inputBinding:
      position: 5
      prefix: "--sasamp"

  treads:
    type:
      - int?
    inputBinding:
      position: 5
      prefix: "--threads"


outputs:

  output_dir:
    type:
      - Directory
    outputBinding:
      glob: '*'

baseCommand: [salmon, index]
