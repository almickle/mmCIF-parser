
export default function parse_mmCIF(text) {

    const ATOMS = text.split('\n').filter((line) => line.startsWith('ATOM')).map(array => array.split(' ').filter(element => element !== ''))
    const customLabels = ['', 'id', 'atom', 'atom_type', '', 'residue', 'chain', 'entity_index', 'residue_index', '', 'x', 'y', 'z', 'occupancy', 'isotropic_temperature_factor', 'formal_charge', 'author_residue_index', 'author_residue', 'author_chain', 'author_type', '' ]    
    const newAtoms = ATOMS.map((array) => {
        return (
            array.map((entry, index) => {
                const label = customLabels[index]
                if(label !== '')
                return (
                    {[label]: entry}
                )}).filter((entry) => entry)
        )
    })

    return newAtoms
}
