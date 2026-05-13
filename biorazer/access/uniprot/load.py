from biorazer.sequence.protein import ProteinSequence

class UniprotSequence(ProteinSequence):
    
    @classmethod
    def load(cls, uniprot_id: str) -> "UniprotSequence":
        """
        Import a protein sequence from UniProt database using its ID.

        Parameters
        ----------
        uniprot_id : str
            The UniProt ID of the protein sequence to be imported.

        Returns
        -------
        UniprotSequence
            An instance of `UniprotSequence` containing the imported sequence.
        """
        # Here you would implement the logic to fetch the sequence from UniProt
        # using the provided uniprot_id, and then create an instance of UniprotSequence.
        # This is a placeholder implementation.
        
        # Example: Fetching sequence data (this is just a mockup)
        sequence_data = "MKTAYIAKQRQISFVKSHFSRQDILDLWIYHTQGYFPDWQNYTPGPGIRYPLTFG"
        
        return cls(sequence_data)