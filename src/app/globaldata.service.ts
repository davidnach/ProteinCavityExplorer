export class GlobalData{


	pid : string;
        experimentId : string;
	pocketToResidues;
		
	constructor(){
		this.pocketToResidues = new Map<number,number[]>();	
	}


	setPid(pid : string){
		this.pid = pid;
	}

	getPid(){
		return this.pid;
	}

	setExperimentId(experimentId : string){
		this.experimentId = experimentId;
	}

	getExperimentId(){
		return this.experimentId;

	}

	getPocketResidues(pocketNum : number){
		return this.pocketToResidues[pocketNum];
	}
	
	setPockets(pockets){
		var residuesString;
		var residues = [];
		var i;
		var j;
		for(i = 0; i < pockets.length; i++){
			console.log(pockets[i].residues);
			residuesString = pockets[i].residues.split(',');
			for (j = 0; j < residuesString.length; j++){
				residues.push(Number(residuesString[j]));
			}
			this.pocketToResidues.set(i,residues);
			residues = [];
		}
	}
}







