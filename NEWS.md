- v0.2.0
  - Add in parsenewick() to parse newick trees using Tokenize
  - Rename iterator and filter functions to nodeiter(),
     branchnamefilter(), etc and use function first ordering to allow
     do block.
  - Major bugfixes and addition of extensive testing of code.
- v0.1.3
  - show() now working for trees again
  - Now provides implementation of getparent(), getancestors(), getchildren(), getdescendants()
  - update docs
- v0.1.2
  - @rput and @rget now work with new RCall RClass interface
  - Can generate random ultrametric trees (using Ultrametric)
  - Can translate and generate BinaryTrees of any kind where data types have default constructors, not just NamedTrees
  - various bugfixes in interfaces
- v0.1.1
  - R <-> Julia interface
  - OrderedDict means tips remain in order originally presented
- v0.1.0
  - Initial release
