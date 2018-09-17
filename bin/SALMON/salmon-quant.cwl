cwlVersion: v1.0
class: CommandLineTool

requirements:
- class: InlineJavascriptRequirement


hints:
- class: DockerRequirement
  dockerPull: tacazares/sc


inputs:

  libtype:
    type:
     - string
    inputBinding:
      position: 1
      prefix: "-l"

  index_file:
    type:
      - Directory
    inputBinding:
      position: 2
      prefix: "-i"

  fastq1_file:
    type:
      - File
    inputBinding:
      position: 6
      prefix: "--mates1"

  fastq2_file:
    type:
      - File
    inputBinding:
      position: 7
      prefix: "--mates2"

  output_filename:
    type:
      - string
    inputBinding:
      position: 8
      prefix: "-o"

  sequence_bias:
    type:
     - boolean
    inputBinding:
      position: 9
      prefix: "--seqBias"

  gc_bias:
    type:
     - boolean
    inputBinding:
      position: 10
      prefix: "--gcBias"

outputs:

  output_dir:
    type:
      - Directory
    outputBinding:
      glob: '*'

baseCommand: [salmon, quant]
