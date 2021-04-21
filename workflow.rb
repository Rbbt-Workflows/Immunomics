require 'rbbt-util'
require 'rbbt/workflow'

Misc.add_libdir if __FILE__ == $0

#require 'rbbt/sources/Immunomics'

module Immunomics
  extend Workflow

end

require 'Immunomics/tasks/epitopes'
require 'Immunomics/tasks/processing'
require 'Immunomics/tasks/presentation'

#require 'rbbt/knowledge_base/Immunomics'
#require 'rbbt/entity/Immunomics'

if defined? Sample
  require 'Immunomics/tasks/sample'
end
