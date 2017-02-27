<TeXmacs|1.99.5>

<style|generic>

<\body>
  <doc-data|<doc-title|On the implementation of symmetric
  extensions>|<doc-author|<author-data|<author-name|Denis
  Rosset>>>|<doc-date|February 23, 2017>>

  <paragraph*|Notation. \V>We work with finite Hilbert spaces. For a Hilbert
  space <math|\<cal-H\>> of dimension <math|d>,
  <math|<math-ss|Lin><around*|(|\<cal-H\>|)>> is the set of linear operators
  on <math|\<cal-H\>>, while <math|<math-ss|Herm><around*|(|\<cal-H\>|)>=<around*|{|\<sigma\>\<in\><math-ss|Lin><around*|(|\<cal-H\>|)><text|
  s.t. >\<sigma\>=\<sigma\><rsup|\<ast\>>|}>> is the set of Hermitian
  operators, where <math|\<sigma\><rsup|\<ast\>>> denotes the adjoint of
  <math|\<sigma\>>. Semidefinite positive operators are written
  <math|<math-ss|Herm><rsub|+><around*|(|\<cal-H\>|)>=<around*|{|\<sigma\>\<in\><math-ss|Herm><around*|(|\<cal-H\>|)><text|
  s.t. <math|\<forall\><around*|\||\<varphi\>|\<rangle\>>\<in\>\<cal-H\>,>><around*|\<langle\>|\<varphi\><mid|\|>\<sigma\><mid|\|>\<varphi\>|\<rangle\>>\<geqslant\>0|}>>.

  This last condition, <math|<text|<math|\<forall\><around*|\||\<varphi\>|\<rangle\>>\<in\>\<cal-H\>,>><around*|\<langle\>|\<varphi\><mid|\|>\<sigma\><mid|\|>\<varphi\>|\<rangle\>>\<geqslant\>0>
  is written more compactly <math|\<sigma\>\<succcurlyeq\>0>.

  By the existence of the finite computational basis, we can write
  indifferently <math|<around*|\||\<varphi\>|\<rangle\>>\<in\>\<cal-H\>> or
  <math|<wide|\<varphi\>|\<vect\>>\<in\><with|math-font|Bbb*|C><rsup|d>>; we
  also write <math|X\<in\><math-ss|Lin><around*|(|\<cal-H\>|)>> or
  <math|X\<in\><with|math-font|Bbb*|C><rsup|d\<times\>d>>.

  <section|Approximations of the separable cone>

  Let <math|\<cal-H\><rsub|<text|A>>=<with|math-font|Bbb*|C><rsup|d<rsub|<text|A>>>>
  and <math|\<cal-H\><rsub|<text|B>>=<with|math-font|Bbb*|C><rsup|d<rsub|<text|B>>>>
  two finite Hilbert spaces of dimension <math|d<rsub|<text|A>>>,
  <math|d<rsub|<text|B>>>.

  <\definition>
    The <em|separable cone> <math|<math-ss|Sep><around*|(|\<cal-H\><rsub|<text|A>><mid|\|>\<cal-H\><rsub|<text|B>>|)>>
    contains the Hermitian operators <math|\<sigma\>\<in\><math-ss|Herm><around*|(|\<cal-H\><rsub|<text|A>>\<otimes\>\<cal-H\><rsub|<text|B>>|)>>
    that possess a separable decomposition: there is a finite number of
    weights <math|p<rsub|i>\<geqslant\>0>, states
    <math|<around*|\||\<alpha\><rsub|i>|\<rangle\>>\<in\>\<cal-H\><rsub|<text|A>>>
    and <math|<around*|\||\<beta\><rsub|i>|\<rangle\>>\<in\>\<cal-H\><rsub|<text|B>>>
    such that

    <\equation>
      \<sigma\>=<big|sum><rsub|i>p<rsub|i><around*|\||\<alpha\><rsub|i>|\<rangle\>><around*|\<langle\>|\<alpha\><rsub|i>|\|>\<otimes\><around*|\||\<beta\><rsub|i>|\<rangle\>><around*|\<langle\>|\<beta\><rsub|i>|\|>.
    </equation>
  </definition>

  Separable states, for example, are such that
  <math|\<rho\>\<in\><math-ss|Sep><around*|(|\<cal-H\><rsub|<text|A>><mid|\|>\<cal-H\><rsub|<text|B>>|)>>
  and <math|tr<around*|[|\<rho\>|]>=1>; however, our definition of
  <math|<math-ss|Sep><around*|(|\<cal-H\><rsub|<text|A>><mid|\|>\<cal-H\><rsub|<text|B>>|)>>
  also include subnormalized states, or states with
  <math|tr<around*|[|\<rho\>|]>\<gtr\>1>.

  An entanglement witness is an operator <math|W\<in\><math-ss|Herm><around*|(|\<cal-H\><rsub|<text|A>>\<otimes\>\<cal-H\><rsub|<text|B>>|)>>
  such that <math|tr<around*|[|W*\<sigma\>|]>\<geqslant\>0> for all
  <math|\<sigma\>\<in\><math-ss|Sep><around*|(|\<cal-H\><rsub|<text|A>><mid|\|>\<cal-H\><rsub|<text|B>>|)>>.
  Thus, the cone of all such operators is the dual to
  <math|<math-ss|Sep><around*|(|\<cal-H\><rsub|<text|A>><mid|\|>\<cal-H\><rsub|<text|B>>|)>>.

  However, deciding whether \P<math|\<rho\>\<in\><math-ss|Sep><around*|(|\<cal-H\><rsub|<text|A>><mid|\|>\<cal-H\><rsub|<text|B>>|)>>\Q
  is hard<nbsp><cite|Gurvits2004>. Instead, we consider an outer
  approximation <math|<wide|<math-ss|S>|~>\<supseteq\><math-ss|Sep><around*|(|\<cal-H\><rsub|<text|A>><mid|\|>\<cal-H\><rsub|<text|B>>|)>>
  containing all <math|\<rho\>> possessing a <math|k>-symmetric extension
  <math|\<tau\>>, where <math|k> is the number of copies of the
  <math|\<cal-H\><rsub|<text|B>>> considered (precise definition below).

  Optionally, we can require <math|\<tau\>> to have positive partial
  transpose (PPT), where the transpose is taken over a number <math|c> of
  copies of the B subsystem (a precise definition will follow). Several PPT
  constraints can be present at the same time, each corresponding to a number
  <math|1\<leqslant\>c<rsub|i>\<leqslant\>k>.

  The approximation <math|<wide|<math-ss|S>|~>> is fully specified by the
  number <math|k> of copies of <math|\<cal-H\><rsub|<text|B>>> and the set
  <math|C=<around*|{|c<rsub|i>|}>> defining the PPT constraints.

  <paragraph*|Definitions in the litterature. \V>The original paper by
  Doherty et al.<nbsp><cite|Doherty2004> uses either <math|C=\<emptyset\>> or
  the full set of PPT constraints <math|C=<around*|{|1,\<ldots\>,k|}>>.
  Subsequent works by Navascués et al.<nbsp><cite|Navascues2009|Navascues2009a>,
  however, use either <math|C=\<emptyset\>> or
  <math|C=<around*|{|ceil<around*|(|k/2|)>|}>> where
  <math|ceil<around*|(|x|)>> is the smallest integer <math|\<geqslant\>x>.
  <strong|VERIFY> This incomplete PPT constraint is motivated by the
  definition of <em|inner> approximations of
  <math|<math-ss|Sep><around*|(|\<cal-H\><rsub|<text|A>><mid|\|>\<cal-H\><rsub|<text|B>>|)>>,
  which we do not study in the present note.

  The toolbox QETLAB<nbsp><cite|Johnston2016> uses this last convention, with
  either <math|C=\<emptyset\>> or <math|C=<around*|{|ceil<around*|(|k/2|)>|}>>.
  <strong|VERIFY>

  We now precise the required notions.

  <subsection|Mathematical definitions>

  <\definition>
    <dueto|Doherty<nbsp><cite|Doherty2004>><label|Def:SymmetricExtension>The
    operator <math|\<sigma\>\<in\><math-ss|Herm><around*|(|\<cal-H\><rsub|<text|A>>\<otimes\>\<cal-H\><rsub|<text|B>>|)>>
    has a <math|k>-symmetric extension if there exists a semidefinite
    positive operator <math|\<tau\>\<in\><math-ss|Herm><rsub|+><around*|(|\<cal-H\><rsub|<text|A>>\<otimes\>\<cal-H\><rsub|<text|B><rsub|1>>\<otimes\>\<ldots\>\<otimes\>\<cal-H\><rsub|<text|B><rsub|k>>|)>>
    such that:

    <\itemize-minus>
      <item>when restricted to <math|\<cal-H\><rsub|<text|A>>\<otimes\>\<cal-H\><rsub|<text|B><rsub|1>>>
      by the partial trace, <math|\<tau\>> reproduces <math|\<sigma\>>:

      <\equation>
        tr<rsub|<text|B><rsub|2>\<ldots\><text|B><rsub|k>><around*|[|\<tau\>|]>=\<sigma\>,
      </equation>

      and the partial trace is defined by

      <\equation>
        tr<rsub|<text|B><rsub|2>\<ldots\><text|B><rsub|k>><around*|[|\<tau\>|]>=<big|sum><rsub|j<rsub|2>\<ldots\>j<rsub|k>><around*|[|<with|math-font|Bbb*|1><rsub|<text|A>>\<otimes\><with|math-font|Bbb*|1><rsub|<text|B><rsub|1>>\<otimes\><around*|\<langle\>|j<rsub|2>\<ldots\>j<rsub|k>|\|>|]>*\<tau\>
        <around*|[|<with|math-font|Bbb*|1><rsub|<text|A>>\<otimes\><with|math-font|Bbb*|1><rsub|<text|B><rsub|1>>\<otimes\><around*|\||j<rsub|2>\<ldots\>j<rsub|k>|\<rangle\>>|]>,
      </equation>

      <item>symmetry, first variant: <math|\<tau\>> is symmetric under all
      permutations of <math|\<cal-H\><rsub|<text|B><rsub|1>>\<otimes\>\<ldots\>\<otimes\>\<cal-H\><rsub|<text|B><rsub|k>>>,

      <item>symmetry, second variant: the support and range of <math|\<tau\>>
      is entirely contained in <math|\<cal-H\><rsub|<text|A>>> and the
      symmetric subspace of <math|\<cal-H\><rsub|<text|B><rsub|1>>\<otimes\>\<ldots\>\<otimes\>\<cal-H\><rsub|<text|B><rsub|k>>>.
    </itemize-minus>
  </definition>

  Let us clarify the notions used in this definition. The first variant uses
  the following.

  <\definition>
    Let <math|\<pi\>\<in\>\<cal-S\><rsub|k>> be an element of the symmetric
    group of degree <math|k>. Alternatively, let
    <math|<around*|{|\<pi\><rsub|1>,\<ldots\>,\<pi\><rsub|k>|}>> be a set
    such that the integers <math|1.\<ldots\>k> appear exactly once.

    The permutation operator <math|P<rsub|\<pi\>>> associated to
    <math|\<pi\>> is given by:

    <\equation>
      P<rsub|\<pi\>>=<big|sum><rsub|i j<rsub|1>\<ldots\>j<rsub|k>><around*|\||i
      j<rsub|\<pi\><rsub|1>>j<rsub|\<pi\><rsub|2>>\<ldots\>j<rsub|\<pi\><rsub|k>>|\<rangle\>><around*|\<langle\>|i
      j<rsub|1>j<rsub|2>\<ldots\>j<rsub|k>|\|>,
    </equation>

    where the sum is done over <math|i\<in\><around*|{|1,\<ldots\>,d<rsub|<text|A>>|}>>
    and <math|j<rsub|1>,\<ldots\>,j<rsub|k>\<in\><around*|{|1,\<ldots\>,d<rsub|<text|B>>|}>>.
  </definition>

  Then <math|\<tau\>> is symmetric under all permutations of
  <math|\<cal-H\><rsub|<text|B><rsub|1>>\<otimes\>\<ldots\>\<otimes\>\<cal-H\><rsub|<text|B><rsub|k>>>
  if <math|P<rsub|\<pi\>>*\<tau\>=\<tau\>*P<rsub|\<pi\>>> for all
  <math|\<pi\>\<in\>\<cal-S\><rsub|k>>.

  The second variant uses the following.

  <\definition>
    <label|Def:SymmetricSubspace>The symmetric subspace of
    <math|\<cal-H\><rsub|<text|B><rsub|1>>\<otimes\>\<ldots\>\<otimes\>\<cal-H\><rsub|<text|B><rsub|k>>>
    has a basis composed of elements of the form:

    <\equation>
      <frac|1|cte><big|sum><rsub|\<pi\>\<in\>\<cal-S\><rsub|k>><around*|\||j<rsub|\<pi\><rsub|1>>j<rsub|\<pi\><rsub|2>>\<ldots\>j<rsub|\<pi\><rsub|k>>|\<rangle\>>,<space|2em>1\<leqslant\>j<rsub|1>\<leqslant\>j<rsub|2>\<leqslant\>\<ldots\>\<leqslant\>j<rsub|k>\<leqslant\>d<rsub|<text|B>>.
    </equation>

    We write these elements <math|<around*|\||\<omega\><rsub|\<ell\>>|\<rangle\>>>,
    with the index <math|\<ell\>=1,\<ldots\>,L> running over all
    <math|j<rsub|1>\<ldots\>j<rsub|k>> satisfying the order condition above.
    For each <math|<around*|\||\<omega\><rsub|\<ell\>>|\<rangle\>>>, the
    constant is chosen such that the resulting coefficients are 0/1.
  </definition>

  <\proposition>
    An operator <math|\<tau\>\<in\><math-ss|Lin><around*|(|\<cal-H\><rsub|A>\<otimes\>\<cal-H\><rsub|<text|B><rsub|1>>\<otimes\>\<ldots\>\<otimes\>\<cal-H\><rsub|<text|B><rsub|k>>|)>>
    obeys the \Psupport and range\Q condition of
    Definition<nbsp><reference|Def:SymmetricExtension> if <math|\<tau\>> can
    be written:

    <\equation>
      \<tau\>=<big|sum><rsub|i<rsub|1>i<rsub|2>\<ell\><rsub|1>\<ell\><rsub|2>>\<beta\><rsub|i<rsub|1>i<rsub|2>\<ell\><rsub|1>\<ell\><rsub|2>><around*|\||i<rsub|1>|\<rangle\>><around*|\<langle\>|i<rsub|2>|\|>\<otimes\><around*|\||\<omega\><rsub|\<ell\><rsub|1>>|\<rangle\>><around*|\<langle\>|\<omega\><rsub|\<ell\><rsub|2>>|\|>
    </equation>

    for complex coefficients <math|\<beta\><rsub|i<rsub|1>i<rsub|2>\<ell\><rsub|1>\<ell\><rsub|2>>>.
  </proposition>

  <subsection|PPT constraints>

  For <math|1\<leqslant\>c\<leqslant\>k>, we define the partial transpose
  operation <math|\<top\><rsub|c>> such that
  <math|\<top\><rsub|c><around*|[|\<tau\>|]>=\<tau\><rsup|\<top\><rsub|<text|B><rsub|1>>\<ldots\>\<top\><rsub|<text|B><rsub|c>>>>.
  Precisely, if

  <\equation>
    \<tau\>=<big|sum><rsub|i i<rprime|'>j<rsub|1>j<rsub|1><rprime|'>\<ldots\>j<rsub|k>j<rsub|k><rprime|'>>\<beta\><rsub|i
    i<rprime|'>j<rsub|1>j<rsub|1><rprime|'>\<ldots\>j<rsub|k>j<rsub|k><rprime|'>><around*|\||i
    j<rsub|1>\<ldots\>j<rsub|k>|\<rangle\>><around*|\<langle\>|i<rprime|'>
    j<rprime|'><rsub|1>\<ldots\>j<rsub|k><rprime|'>|\|>,
  </equation>

  then

  <\equation>
    \<top\><rsub|c><around*|[|\<tau\>|]>=<big|sum><rsub|i
    i<rprime|'>j<rsub|1>j<rsub|1><rprime|'>\<ldots\>j<rsub|k>j<rsub|k><rprime|'>>\<beta\><rsub|i
    i<rprime|'>j<rsub|1>j<rsub|1><rprime|'>\<ldots\>j<rsub|k>j<rsub|k><rprime|'>><around*|\||i
    j<rprime|'><rsub|1>\<ldots\>j<rprime|'><rsub|c>j<rsub|c+1>\<ldots\>j<rsub|k>|\<rangle\>><around*|\<langle\>|i<rprime|'>
    j<rsub|1>\<ldots\>j<rsub|c>j<rsub|c+1><rprime|'>\<ldots\>j<rsub|k><rprime|'>|\|>.
  </equation>

  The positive partial transpose constraints are given by
  <math|\<top\><rsub|c><around*|[|\<tau\>|]>\<succcurlyeq\>0> for each
  <math|c\<in\>C>.

  <section|Semidefinite formulation using equality constraints>

  This form is called the primal form in the documentation of most
  primal-dual solvers such as SeDuMi<nbsp><cite|Sturm2002>,
  SDPT3<nbsp><cite|Tutuncu2003>, etc...

  It is <em|not> the default form provided by YALMIP, however, the setting
  <verbatim|dualize> can be set to <math|1> to employ this form. It
  <em|seems> to be the form preferred by CVX version 2 (<strong|verify>);
  however, CVX version 3 can formulate problems in either the primal or the
  dual form; this is controlled by the command <verbatim|cvx_dualize>.

  <subsection|Without using the \Psupport and range\Q condition>

  <paragraph*|Primal variables. \V>Our variables are
  <math|\<sigma\>\<in\><math-ss|Herm><rsub|+><around*|(|\<cal-H\><rsub|<text|A>>\<otimes\>\<cal-H\><rsub|<text|B>>|)>>,
  the element of our approximation <math|<wide|<math-ss|S>|~>>, and the
  extension <math|\<tau\>\<in\><math-ss|Herm><rsub|+><around*|(|\<cal-H\><rsub|<text|A>>\<otimes\>\<cal-H\><rsub|<text|B><rsub|1>>\<otimes\>\<ldots\>\<otimes\>\<cal-H\><rsub|<text|B><rsub|k>>|)>>.

  <paragraph*|Equality constraints. \V>The extension reproduces the bipartite
  state:<math|<wide|asdasda|^>>

  <\equation>
    tr<rsub|<text|B><rsub|2>\<ldots\><text|B><rsub|k>><around*|[|\<tau\>|]>=\<sigma\>,
  </equation>

  and the extension is symmetric:

  <\equation>
    <label|Eq:InSymmetricSubspace>\<Pi\>\<tau\>\<Pi\>=\<tau\>,
  </equation>

  where <math|\<Pi\>=<big|sum><rsub|\<pi\>>P<rsub|\<pi\>>> is the projection
  onto the symmetric subspace.

  <paragraph*|Remarks. \V>The semidefinite constraints
  <math|\<tau\>\<succcurlyeq\>0> and <math|\<sigma\>\<succcurlyeq\>0> are
  automatically present in the primal form and do not need to be added
  manually.\ 

  When the symmetric extension is contained in a larger semidefinite program,
  the variable <math|\<sigma\>> can well be represented by linear
  combinations of other variables. In that case, it is not necessary to
  define an additional semidefinite variable <math|\<sigma\>>; instead, the
  equality constraints can be modified accordingly.

  <subsubsection|PPT constraints>

  <paragraph*|Additional primal variables. \V>Each PPT constraint
  <math|c\<in\>C> requires a new variable
  <math|\<rho\><rsub|c>\<in\><math-ss|Herm><rsub|+><around*|(|\<cal-H\><rsub|<text|A>>\<otimes\>\<cal-H\><rsub|<text|B><rsub|1>>\<otimes\>\<ldots\>\<otimes\>\<cal-H\><rsub|<text|B><rsub|k>>|)>>.

  <paragraph*|Additional equality constraints. \V>We require:

  <\equation>
    \<rho\><rsub|c>=\<top\><rsub|c><around*|[|\<tau\>|]>.
  </equation>

  <subsubsection|Final remarks>

  While conceptually simple, this formulation is not particularly efficient
  due to the high number of equality constraints
  in<nbsp><eqref|Eq:InSymmetricSubspace> and the large dimension of the
  matrices. The next section provides a more efficient formulation.

  <subsection|With the \Psupport and range\Q condition>

  We use the following parameterization of <math|\<tau\>>, when its support
  and range is restricted to the symmetric subspace.

  <\proposition>
    <label|Prop:SupportAndRange>Let <math|L> be the size of the symmetric
    subspace in Definition<nbsp><reference|Def:SymmetricSubspace>, and we
    write <math|\<cal-H\><rsub|<text|L>>=<with|math-font|Bbb*|C><rsup|L>>.
    Then <math|\<tau\>> obeys the \Psupport and range\Q condition if there
    exists <math|<wide|\<tau\>|\<invbreve\>>\<in\><math-ss|Herm><around*|(|\<cal-H\><rsub|<text|A>>\<otimes\>\<cal-H\><rsub|<text|L>>|)>>
    such that

    <\equation>
      \<tau\>=<around*|(|<with|math-font|Bbb*|1><rsub|<text|A>>\<otimes\>\<Omega\>|)>*<wide|\<tau\>|\<invbreve\>>*<around*|(|<with|math-font|Bbb*|1><rsub|<text|A>>\<otimes\>\<Omega\>|)><rsup|\<ast\>>,
    </equation>

    where <math|\<Omega\>> represents the change of basis from
    <math|\<cal-H\><rsub|<text|L>>> to <math|\<cal-H\><rsub|<text|B><rsub|1>>\<otimes\>\<ldots\>\<otimes\>\<cal-H\><rsub|<text|B><rsub|k>>>:

    <\equation>
      \<Omega\>=<big|sum><rsub|\<ell\>><around*|\||\<omega\><rsub|\<ell\>>|\<rangle\>><around*|\<langle\>|\<ell\>|\|>\<in\><with|math-font|Bbb*|C><rsup|L\<times\>d<rsub|<text|B>><rsup|k>>.
    </equation>

    Note: <math|\<Omega\>> is obtained by simply concatening the column
    vectors corresponding to <math|<around*|\||\<omega\><rsub|\<ell\>>|\<rangle\>>>.
  </proposition>

  \;

  <paragraph*|Primal variables. \V>Our primal variables are
  <math|\<sigma\>\<in\><math-ss|Herm><rsub|+><around*|(|\<cal-H\><rsub|<text|A>>\<otimes\>\<cal-H\><rsub|<text|B>>|)>>,
  the element of our approximation <math|<wide|<math-ss|S>|~>>, and the
  parameterization <math|<wide|\<tau\>|\<invbreve\>>\<in\><math-ss|Herm><rsub|+><around*|(|\<cal-H\><rsub|<text|A>>\<otimes\>\<cal-H\><rsub|<text|L>>|)>>.
  Positive semidefinitiveness of both is implied by the primal formulation.
  Note that <math|<wide|\<tau\>|\<invbreve\>>\<succcurlyeq\>0> implies
  <math|\<tau\>\<succcurlyeq\>0>.

  <paragraph*|Equality constraints. \V>We have a single equality constraint,
  that the extension reproduces the bipartite state:

  <\equation>
    tr<rsub|<text|B><rsub|2>\<ldots\><text|B><rsub|k>><around*|[|<around*|(|<with|math-font|Bbb*|1><rsub|<text|A>>\<otimes\>\<Omega\>|)>*<wide|\<tau\>|\<invbreve\>>*<around*|(|<with|math-font|Bbb*|1><rsub|<text|A>>\<otimes\>\<Omega\>|)>|]>=\<sigma\>.
  </equation>

  <subsubsection|PPT constraints>

  Each PPT constraint <math|c\<in\>C> breaks the symmetry of
  <math|\<cal-H\><rsub|<text|B><rsub|1>>\<otimes\>\<ldots\>\<otimes\>\<cal-H\><rsub|<text|B><rsub|k>
  >>; instead, <math|\<rho\><rsub|c>> has support and range in the symmetric
  subspaces of <math|><math|\<cal-H\><rsub|<text|B><rsub|1>>\<otimes\>\<ldots\>\<otimes\>\<cal-H\><rsub|<text|B><rsub|c>>>
  and <math|\<cal-H\><rsub|<text|B><rsub|c+1>>\<otimes\>\<ldots\>\<otimes\>\<cal-H\><rsub|<text|B><rsub|k>>>.
  We reuse the trick of Proposition<nbsp><reference|Prop:SupportAndRange> by
  introducing the bases <math|<around*|{|<around*|\||\<omega\><rsub|\<ell\>><rprime|'>|\<rangle\>>|}>>
  and <math|<around*|{|<around*|\||\<omega\><rsub|\<ell\>><rprime|''>|\<rangle\>>|}>>
  of these respective symmetric subspaces, with respective dimension
  <math|L<rprime|'>> and <math|L<rprime|''>>.

  <paragraph*|Additional primal variables. \V>The additional variable
  <math|<wide|\<rho\>|\<invbreve\>><rsub|c>\<in\><math-ss|Herm><rsub|+><around*|(|\<cal-H\><rsub|<text|A>>\<otimes\>\<cal-H\><rsub|<text|L><rprime|'>>\<otimes\>\<cal-H\><rsub|<text|L><rprime|''>>|)>>
  parameterizes

  <\equation>
    \<rho\><rsub|c>=<around*|(|<with|math-font|Bbb*|1><rsub|<text|A>>\<otimes\>\<Omega\><rprime|'>\<otimes\>\<Omega\><rprime|''>|)>*<wide|\<rho\>|\<invbreve\>><rsub|c>*<around*|(|<with|math-font|Bbb*|1><rsub|<text|A>>\<otimes\>\<Omega\><rprime|'>\<otimes\>\<Omega\><rprime|''>|)><rsup|\<ast\>>,
  </equation>

  where <math|\<Omega\><rprime|'>> represents the change of basis from
  <math|\<cal-H\><rsub|<text|L><rprime|'>>> to
  <math|\<cal-H\><rsub|<text|B><rsub|1>>\<otimes\>\<ldots\>\<otimes\>\<cal-H\><rsub|<text|B><rsub|c>>>:

  <\equation>
    \<Omega\><rprime|'>=<big|sum><rsub|\<ell\><rprime|'>><around*|\||\<omega\><rprime|'><rsub|\<ell\><rprime|'>>|\<rangle\>><around*|\<langle\>|\<ell\><rprime|'>|\|>\<in\><with|math-font|Bbb*|C><rsup|L<rprime|'>\<times\>d<rsub|<text|B>><rsup|c>>,
  </equation>

  and <math|\<Omega\><rprime|''>> represents the change of basis from
  <math|\<cal-H\><rsub|<text|L><rprime|''>>> to
  <math|\<cal-H\><rsub|<text|B><rsub|c+1>>\<otimes\>\<ldots\>\<otimes\>\<cal-H\><rsub|<text|B><rsub|k>>>:

  <\equation>
    \<Omega\><rprime|''>=<big|sum><rsub|\<ell\><rprime|''>><around*|\||\<omega\><rprime|''><rsub|\<ell\><rprime|''>>|\<rangle\>><around*|\<langle\>|\<ell\><rprime|''>|\|>\<in\><with|math-font|Bbb*|C><rsup|L<rprime|''>\<times\>d<rsub|<text|B>><rsup|k-c>>.
  </equation>

  <paragraph*|Additional equality constraints. \V>We require:

  <\equation>
    \<rho\><rsub|c>=\<top\><rsub|c><around*|[|<around*|(|<with|math-font|Bbb*|1><rsub|<text|A>>\<otimes\>\<Omega\><rprime|'>\<otimes\>\<Omega\><rprime|''>|)>*<wide|\<rho\>|\<invbreve\>><rsub|c>*<around*|(|<with|math-font|Bbb*|1><rsub|<text|A>>\<otimes\>\<Omega\><rprime|'>\<otimes\>\<Omega\><rprime|''>|)><rsup|\<ast\>>|]>.
  </equation>

  Note: implemented naively, this adds many redudant equality constraints,
  causing numerical instability. Instead, write, for
  <math|i<rsub|1>,i<rsub|2>=1,\<ldots\>,d<rsub|<text|A>>> and
  <math|\<ell\><rsub|1><rprime|'>,\<ell\><rsub|2><rprime|'>=1,\<ldots\>,L<rprime|'>>
  and <math|\<ell\><rprime|''><rsub|1>,\<ell\><rsub|2><rprime|''>=1,\<ldots\>,L<rprime|''>>:

  <\equation>
    <around*|\<langle\>|i<rsub|1>\<ell\><rsub|2><rprime|'>\<ell\><rsub|1><rprime|''><mid|\|><wide|\<rho\>|\<invbreve\>><rsub|c><mid|\|>i<rsub|2>\<ell\><rsub|1><rprime|'>\<ell\><rsub|2><rprime|''>|\<rangle\>>=<around*|\<langle\>|i<rsub|1>\<ell\><rsub|1><mid|\|><wide|\<tau\>|\<invbreve\>><mid|\|>i<rsub|2>\<ell\><rsub|2>|\<rangle\>>.
  </equation>

  We left undefined the relations <math|\<ell\><rsub|1>=f<around*|(|\<ell\><rprime|'><rsub|2>,\<ell\><rsub|1><rprime|''>|)>>
  and <math|\<ell\><rsub|2>=f<around*|(|\<ell\><rsub|1><rprime|'>,\<ell\><rsub|2><rprime|''>|)>>.

  The function <math|f> is defined as <math|\<ell\>\<equiv\>f<around*|(|\<ell\><rprime|'>,\<ell\><rprime|''>|)>>,
  with <math|1\<leqslant\>\<ell\>\<leqslant\>L>,
  <math|1\<leqslant\>\<ell\><rprime|'>\<leqslant\>L<rprime|'>>,
  <math|1\<leqslant\>\<ell\><rprime|''>\<leqslant\>L<rprime|''>>. To obtain
  <math|\<ell\>>, we first decompose take a nonzero
  <math|<around*|\||j<rsub|1>\<ldots\>j<rsub|c>|\<rangle\>>> present in
  <math|<around*|\||\<omega\><rprime|'><rsub|\<ell\><rprime|'>>|\<rangle\>>>,
  and a nonzero <math|<around*|\||j<rsub|c+1>\<ldots\>j<rsub|k>|\<rangle\>>>
  present in <math|<around*|\||\<omega\><rsub|\<ell\><rprime|''>><rprime|''>|\<rangle\>>>.
  The basis <math|<around*|\||\<omega\><rsub|\<ell\>>|\<rangle\>>> containing
  a nonzero <math|<around*|\||j<rsub|1>\<ldots\>j<rsub|k>|\<rangle\>>> gives
  the required <math|\<ell\>>.

  <section|Semidefinite formulation using inequality constraints>

  Let <math|<around*|{|\<alpha\><rsub|i>|}>\<in\><math-ss|Herm><around*|(|\<cal-H\><rsub|<text|A>>|)>>
  and <math|<around*|{|\<beta\><rsub|j>|}>\<in\><math-ss|Herm><around*|(|\<cal-H\><rsub|<text|B>>|)>>
  be bases of their respective spaces, such that
  <math|\<alpha\><rsub|1>=<with|math-font|Bbb*|1><rsub|<text|A>>> and
  <math|\<beta\><rsub|1>=<with|math-font|Bbb*|1><rsub|<text|B>>>, and
  <math|i=1,\<ldots\>,d<rsub|<text|A><rsup|>><rsup|2>> while
  <math|j=1,\<ldots\>,d<rsub|<text|B>><rsup|2>>. These bases should be
  orthogonal, but not necessarily normalized.

  Then, the operator <math|\<sigma\>\<in\><math-ss|Herm><around*|(|\<cal-H\><rsub|<text|A>>\<otimes\>\<cal-H\><rsub|<text|B>>|)>>
  has a decomposition:

  <\equation>
    \<sigma\>=<big|sum><rsub|i j>s<rsub|i
    j><around*|(|\<alpha\><rsub|i>\<otimes\>\<beta\><rsub|j>|)>,
  </equation>

  while its symmetric extension <math|\<tau\>> has the decomposition:

  <\equation>
    \<tau\>=<big|sum><rsub|i j<rsub|1>\<ldots\>j<rsub|k>>t<rsub|i
    j<rsub|1>\<ldots\>j<rsub|k>><around*|(|\<alpha\><rsub|i>\<otimes\>\<beta\><rsub|j<rsub|1>>\<otimes\>\<ldots\>\<otimes\>\<beta\><rsub|j<rsub|k>>|)>,
  </equation>

  with the constraint that <math|t<rsub|i j 1\<ldots\>1>=s<rsub|i j>> and
  that <math|t<rsub|i j<rsub|\<pi\><rsub|1>>\<ldots\>j<rsub|\<pi\><rsub|k>>>=t<rsub|i
  j<rsub|1>\<ldots\>j<rsub|k>>> for all permutations
  <math|\<pi\>\<in\>\<cal-S\><rsub|k>>.

  <paragraph*|Dual variables. \V>The coefficients <math|t<rsub|i
  j<rsub|1>\<ldots\>j<rsub|k>>> are variables of our formulation, with the
  following exceptions:

  <\itemize-minus>
    <item>when the <math|j<rsub|1>\<ldots\>j<rsub|k>> are not sorted in
    increasing order, the coefficient <math|t<rsub|i
    j<rsub|1>\<ldots\>j<rsub|k>>> is an alias for <math|t<rsub|i
    j<rsub|1><rprime|'>\<ldots\>j<rprime|'><rsub|k>>> where
    <math|<around*|(|j<rsub|1><rprime|'>,\<ldots\>,j<rsub|k><rprime|'>|)>>
    are the indices <math|<around*|(|j<rsub|1>,\<ldots\>,j<rsub|k>|)>> sorted
    in increasing order (note: \Pincreasing\Q is actually \Pnondecreasing\Q,
    due to possible repetitions in the index values),

    <item>the coefficient <math|t<rsub|i j<rsub|1>1\<ldots\>1>> is directly
    set to the value <math|s<rsub|i j<rsub|1>>>.
  </itemize-minus>

  These coefficients are all real but otherwise unrestricted. If the
  symmetric extension is part of a larger semidefinite program, the
  coefficients <math|s<rsub|i j>> can be replaced by the relevant (linear)
  expressions.

  Then the constraint:

  <\equation>
    <label|Eq:DualSymPSD>\<tau\>=<big|sum><rsub|i
    j<rsub|1>\<ldots\>j<rsub|k>>t<rsub|i j<rsub|1>\<ldots\>j<rsub|k>><around*|(|\<alpha\><rsub|i>\<otimes\>\<beta\><rsub|j<rsub|1>>\<otimes\>\<ldots\>\<otimes\>\<beta\><rsub|j<rsub|k>>|)>\<geqslant\>0
  </equation>

  correspond to the canonical form of the solver.

  <paragraph*|PPT constraints. \V>These constraints are easy to incorporate.
  For <math|c\<in\>C>, simply add:

  <\equation>
    <big|sum><rsub|i j<rsub|1>\<ldots\>j<rsub|k>>t<rsub|i
    j<rsub|1>\<ldots\>j<rsub|k>><around*|(|\<alpha\><rsub|i>\<otimes\>\<beta\><rsup|\<top\>><rsub|j<rsub|1>>\<otimes\>\<ldots\>\<otimes\>\<beta\><rsub|j<rsub|c>><rsup|\<top\>>\<otimes\>\<beta\><rsub|j<rsub|c+1>>\<otimes\>\<ldots\>\<otimes\>\<beta\><rsub|j<rsub|k>>|)>\<geqslant\>0.
  </equation>

  <subsection|Restriction to the symmetric subspace>

  By an appropriate restriction to the symmetric subspace, the formulation
  above can be improved. The constraint<nbsp><eqref|Eq:DualSymPSD> can be
  imposed only on the subspace defined by <math|<around*|\||i
  j<rsub|1>\<ldots\>j<rsub|k>|\<rangle\>>> when the
  <math|j<rsub|1>\<ldots\>j<rsub|k>> are in increasing order. Note that in
  the present subsection, <math|i> and <math|j<rsub|1>\<ldots\>j<rsub|k>> are
  indices of the bra-ket computational basis, not the indices of the
  Hermitian basis.

  This corresponds to dropping rows and columns of the matrices in the dual
  constraint<nbsp><eqref|Eq:DualSymPSD>. For the PPT constraints, the story
  is similar, except that we require only that the indices
  <math|j<rsub|1>\<ldots\>j<rsub|c>> on the one hand,
  <math|j<rsub|c+1>\<ldots\>j<rsub|k>> on the other hand to be increasing.

  However, this formulation introduce linear dependences in blocks of the
  problem, and can lead to numerical instabilities.

  <\bibliography|bib|tm-plain|~/Documents/Bibliography/ZotOutput.bib>
    <\bib-list|5>
      <bibitem*|1><label|bib-Doherty2004>Andrew<nbsp>C.<nbsp>Doherty,
      Pablo<nbsp>A.<nbsp>Parrilo<localize|, and
      >Federico<nbsp>M.<nbsp>Spedalieri.<newblock> Complete family of
      separability criteria.<newblock> <with|font-shape|italic|Physical
      Review A>, 69(2):22308, feb 2004.<newblock>

      <bibitem*|2><label|bib-Navascues2009>Miguel Navascués, Masaki
      Owari<localize|, and >Martin<nbsp>B.<nbsp>Plenio.<newblock> Complete
      Criterion for Separability Detection.<newblock>
      <with|font-shape|italic|Physical Review Letters>, 103(16):160404, oct
      2009.<newblock>

      <bibitem*|3><label|bib-Navascues2009a>Miguel Navascués, Masaki
      Owari<localize|, and >Martin<nbsp>B.<nbsp>Plenio.<newblock> Power of
      symmetric extensions for entanglement detection.<newblock>
      <with|font-shape|italic|Physical Review A>, 80(5):52306, nov
      2009.<newblock>

      <bibitem*|4><label|bib-Sturm2002>Jos<nbsp>F.<nbsp>Sturm.<newblock>
      Implementation of interior point methods for mixed semidefinite and
      second order cone optimization problems.<newblock>
      <with|font-shape|italic|Optimization Methods and Software>,
      17(6):1105\U1154, jan 2002.<newblock>

      <bibitem*|5><label|bib-Tutuncu2003>R.<nbsp>H.<nbsp>Tütüncü,
      K.<nbsp>C.<nbsp>Toh<localize|, and >M.<nbsp>J.<nbsp>Todd.<newblock>
      Solving semidefinite-quadratic-linear programs using SDPT3.<newblock>
      <with|font-shape|italic|Mathematical Programming>, 95(2):189\U217, feb
      2003.<newblock>
    </bib-list>
  </bibliography>
</body>

<initial|<\collection>
</collection>>

<\references>
  <\collection>
    <associate|Def:SymmetricExtension|<tuple|2|2>>
    <associate|Def:SymmetricSubspace|<tuple|4|2>>
    <associate|Eq:DualSymPSD|<tuple|22|5>>
    <associate|Eq:InSymmetricSubspace|<tuple|10|3>>
    <associate|Prop:SupportAndRange|<tuple|6|4>>
    <associate|auto-1|<tuple|?|1>>
    <associate|auto-10|<tuple|10|3>>
    <associate|auto-11|<tuple|2.1.1|3>>
    <associate|auto-12|<tuple|2.1.1|3>>
    <associate|auto-13|<tuple|2.1.1|3>>
    <associate|auto-14|<tuple|2.1.2|4>>
    <associate|auto-15|<tuple|2.2|4>>
    <associate|auto-16|<tuple|13|4>>
    <associate|auto-17|<tuple|13|4>>
    <associate|auto-18|<tuple|2.2.1|4>>
    <associate|auto-19|<tuple|2.2.1|4>>
    <associate|auto-2|<tuple|1|1>>
    <associate|auto-20|<tuple|17|5>>
    <associate|auto-21|<tuple|3|5>>
    <associate|auto-22|<tuple|21|5>>
    <associate|auto-23|<tuple|22|5>>
    <associate|auto-24|<tuple|3.1|6>>
    <associate|auto-25|<tuple|3.1|?>>
    <associate|auto-3|<tuple|1|1>>
    <associate|auto-4|<tuple|1.1|2>>
    <associate|auto-5|<tuple|1.2|3>>
    <associate|auto-6|<tuple|2|3>>
    <associate|auto-7|<tuple|2.1|3>>
    <associate|auto-8|<tuple|2.1|3>>
    <associate|auto-9|<tuple|2.1|3>>
    <associate|bib-Doherty2004|<tuple|1|?>>
    <associate|bib-Navascues2009|<tuple|2|?>>
    <associate|bib-Navascues2009a|<tuple|3|?>>
    <associate|bib-Sturm2002|<tuple|4|?>>
    <associate|bib-Tutuncu2003|<tuple|5|?>>
  </collection>
</references>

<\auxiliary>
  <\collection>
    <\associate|bib>
      Gurvits2004

      Doherty2004

      Navascues2009

      Navascues2009a

      Johnston2016

      Doherty2004

      Sturm2002

      Tutuncu2003
    </associate>
    <\associate|toc>
      <with|par-left|<quote|4tab>|Notation. \V
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-1><vspace|0.15fn>>

      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|1<space|2spc>Approximations
      of the separable cone> <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-2><vspace|0.5fn>

      <with|par-left|<quote|4tab>|Definitions in the litterature. \V
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-3><vspace|0.15fn>>

      <with|par-left|<quote|1tab>|1.1<space|2spc>Mathematical definitions
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-4>>

      <with|par-left|<quote|1tab>|1.2<space|2spc>PPT constraints
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-5>>

      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|2<space|2spc>Semidefinite
      formulation using equality constraints>
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-6><vspace|0.5fn>

      <with|par-left|<quote|1tab>|2.1<space|2spc>Without using the \Psupport
      and range\Q condition <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-7>>

      <with|par-left|<quote|4tab>|Primal variables. \V
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-8><vspace|0.15fn>>

      <with|par-left|<quote|4tab>|Equality constraints. \V
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-9><vspace|0.15fn>>

      <with|par-left|<quote|4tab>|Remarks. \V
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-10><vspace|0.15fn>>

      <with|par-left|<quote|2tab>|2.1.1<space|2spc>PPT constraints
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-11>>

      <with|par-left|<quote|4tab>|Additional primal variables. \V
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-12><vspace|0.15fn>>

      <with|par-left|<quote|4tab>|Additional equality constraints. \V
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-13><vspace|0.15fn>>

      <with|par-left|<quote|2tab>|2.1.2<space|2spc>Final remarks
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-14>>

      <with|par-left|<quote|1tab>|2.2<space|2spc>With the \Psupport and
      range\Q condition <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-15>>

      <with|par-left|<quote|4tab>|Primal variables. \V
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-16><vspace|0.15fn>>

      <with|par-left|<quote|4tab>|Equality constraints. \V
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-17><vspace|0.15fn>>

      <with|par-left|<quote|2tab>|2.2.1<space|2spc>PPT constraints
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-18>>

      <with|par-left|<quote|4tab>|Additional primal variables. \V
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-19><vspace|0.15fn>>

      <with|par-left|<quote|4tab>|Additional equality constraints. \V
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-20><vspace|0.15fn>>

      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|3<space|2spc>Semidefinite
      formulation using inequality constraints>
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-21><vspace|0.5fn>

      <with|par-left|<quote|4tab>|Dual variables. \V
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-22><vspace|0.15fn>>

      <with|par-left|<quote|4tab>|PPT constraints. \V
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-23><vspace|0.15fn>>

      <with|par-left|<quote|1tab>|3.1<space|2spc>Restriction to the symmetric
      subspace <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-24>>

      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|Bibliography>
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-25><vspace|0.5fn>
    </associate>
  </collection>
</auxiliary>