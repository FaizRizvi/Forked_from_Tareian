cwlVersion: v1.0
class: CommandLineTool

requirements:
- class: InlineJavascriptRequirement

hints:
- class: DockerRequirement
  dockerPull: insilicodb/kallisto

inputs:

  input_filename:
    type:
      - string
    inputBinding:
      position: 5
      prefix: "-i"

  fasta_file:
    type:
      - File

outputs:

  output_dir:
    type:
      - Directory
    outputBinding:
      glob: '*'

baseCommand: [kallisto, index]
