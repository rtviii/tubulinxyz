
Artem, [2025-11-24 16:35]
# Discussed

### Home page layout and behavior:

The common vision seems to be that the landing page should be simple and easily digestable: two general categories of tubulin (dimers and lattices) as well as the universal llm-powered search. 

Here we discussed:
 - possibly collecting queries from users in the longer term for seeing what general motifs people investigate.
 - change top-level classification to just "dimers"(soluble~free a/b) and "lattices" (including full MTs) and not include various oligomers as a category. This will both simplify the initial impression as well as our implicit categorization. 

### Annotations (mutations/modifications/ligand and MAP interactions):

Morisette's dataset of PTMs and Mutations is now fully integrated and browseable via the general tool.
Priority should be "interrogating relationship with disease". Soft focus on human sequences for now.
Add ligand and MAP interactions.
Agreed not to preoccupy ourselves with semantic anotations for now (like associating each PTM/Mutation/Ligand with text and various other databases of disease, celltype, organelles, process etc.; thouhg we can still provide a reference for the user to decide themselves).


### Sequence/structure co-exploration tool:

Pretty much agreed that the current approach is sound enough. It is easy to revisit the master MSAs in the future.

Carsten stresses the need for intuitive colorscheme to highlight only the important features of the loaded data (mutations/ptms, outliers etc.), mentioned possibly supporting fragmented snippets for evolutionary queries?

# Future desirables 

Maxim & Carsten: integrating a recent structure prediction model like Boltz/Collabfold/AF3(unlikely); integrating a structural prediction with an LLM aware of the particulars of tubulin structures through our application's schema would be somewhat of a holy grail in terms of functionality.

# Other

explicitly agreed on some obvious technical choices:
- use [similarity with a] per-family MSA as main classification mechanism for target sequences in the wild
- abandon the idea of any more involved "polymerization state" classification beyond dimers/lattices

# Further steps

- family classification mechanism
- extend the coverage to beta and gamma tubulins
- create a class of ligand, MAP binding site annotations/interactions
- expand filters: integrate sequence queries into the UI
- lots of UI work
