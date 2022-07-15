require 'rbbt-util'
require 'rbbt/workflow'

Misc.add_libdir if __FILE__ == $0

#require 'rbbt/sources/Immunomics'
#
Workflow.require_workflow "Sequence"
Workflow.require_workflow "HTS"
Workflow.require_workflow "Sample"

module Immunomics
  extend Workflow

end

require 'Immunomics/tasks/epitopes'
require 'Immunomics/tasks/processing'
require 'Immunomics/tasks/presentation'
require 'Immunomics/tasks/expression'
require 'Immunomics/tasks/vcf'

#require 'rbbt/knowledge_base/Immunomics'
#require 'rbbt/entity/Immunomics'

if defined? Sample
  require 'Immunomics/tasks/sample'
  require 'Immunomics/tasks/adhoc'
end
