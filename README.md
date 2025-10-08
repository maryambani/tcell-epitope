## T-Cell Epitope Immunogenicity Classifier

This is a personal learning project with the goal to build a binary classifier that can predict whether a T-cell epitope peptide is immunogenic (will provoke a T-cell response from the immune system) or not.

What are T-cells? They are part of the immune system and are able to recognize problematic cells (such as virus-infected cells) by having their receptors bind to a small protein piece (peptide) that is presented on a cell's surface by molecules called the Major Histocompatibility Complex (MHC). The receptors can only bind the peptide-MHC complex and they are specific on what they can and cannot bind to for a particular combination of peptide and MHC. [1]

MHC is highly polymorphic (many different versions exist in nature) giving the immune system the ability to recognize many, many pathogens. Recognizing an epitope and having an MHC-peptide complex is one of the early steps (although not sufficient on its own) to the immune response. [1]

Accurately predicting which antigens will trigger an immune response is important for better understanding autoimmune diseases, to guide vaccine design and to improve drug safety. With the help of ML, predictions can also be used to help researchers identify neo-antigens for personalized cancer immunotherapy. [2]

#### Scope

This will be a supervised binary classification problem. The goal is to train a model on short peptide sequences. For MHC Class I molecules, they tend to fit peptides of about 8-10 amino acids in length, so for sake of simplicity the dataset will start with peptides within that length range, and only the peptide sequence will be considered as the input [3]. For the scope of this project it is acceptable to only consider the peptide sequences, since the peptides labeled as immunogenic from the dataset are inherently experimentally verified in the presence of an MHC molecule, meaning we can implicitly account for MHC binding.

#### Environment setup

`pip install -r requirements.txt`

#### Data sources

The dataset will be created from two complementary datasets: IEDB ([Immune Epitope Database](https://www.iedb.org/result_v3.php?cookie_id=aa2a4e)) for broad T-cell epitopes across different pathogens and self-antigens, and CEDAR ([Cancer Epitope Database and Analysis Resource](https://cedar.iedb.org)) for a subset of specifically cancer-derived epitopes.

#### Data preparation

#### Feature encoding

#### Modeling

#### Evaluation

#### Results

#### References

[1] Cano RLE, Lopera HDE. Introduction to T and B lymphocytes. In: Anaya JM, Shoenfeld Y, Rojas-Villarraga A, et al., editors. Autoimmunity: From Bench to Bedside [Internet]. Bogota (Colombia): El Rosario University Press; 2013 Jul 18. Chapter 5. Available from: https://www.ncbi.nlm.nih.gov/books/NBK459471/

[2] Xu Z, Wang X, Zeng S, Ren X, Yan Y, Gong Z. Applying artificial intelligence for cancer immunotherapy. Acta Pharm Sin B. 2021 Nov;11(11):3393-3405. doi: 10.1016/j.apsb.2021.02.007. Epub 2021 Feb 11. PMID: 34900525; PMCID: PMC8642413.

[3] Alberts B, Johnson A, Lewis J, et al. Molecular Biology of the Cell. 4th edition. New York: Garland Science; 2002. T Cells and MHC Proteins. Available from: https://www.ncbi.nlm.nih.gov/books/NBK26926/
