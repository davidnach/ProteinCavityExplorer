const express = require('express');

const router = express.Router();
const csvToJson = require('csvjson');
const fs = require('fs');
const path = require('path');
var spawn = require("child_process").spawn;
const shell = require('shelljs');




router.post('/createExperiment',  (req,res) => {
    // quick safety check on incoming pid
    if(/^([1-9]{1})([a-zA-Z0-9]{3})$/.test(req.body.pid)) {
	console.log('server printing pid' + req.body.pid);
	var pid = req.body.pid;
	var expID = "";
	var pythonProcess = spawn('python',["pidDirectory.py", pid]);
	pythonProcess.stdout.on('data', function (data){
            // do something with experiment ID (data)
            expID = pid + data;
	    //shell.exec('../../bashrun/runrmutantandvasp.sh ' + pid + ' ' + expID,{shell: '/bin/bash'});
	    console.log('in server');
            var expId = {expId : expID };
	    var data;
	    var pockets;
	    var proteinInfo 
	    try {
	 	data = fs.readFileSync(path.join('/home/completePocketsNRSet',pid + 'A' + '.csv'), {encoding : 'utf8'});
		var options = {
                    delimiter : " ",
                    wrap : true,
                    headers :"pocket area residues"};
		pockets = csvToJson.toObject(data,options);
		console.log(pockets);
		proteinInfo = {expId, pockets};
	    } catch (err) {
		proteinInfo = expId;
		console.log('surf report doesnt exist');
	    }
	    
            res.send(JSON.stringify(proteinInfo));
	});
    } else { res.send('backend validation failed!'); // probably should do better than this
	     // Only malicious users should ever see this, however!
	   }
});



router.post('/retrieveExperiment', (req,res) => {
    console.log(req.body.expId)
    if (/^([1-9]{1})([a-zA-Z0-9]{3})([a-fA-F0-9]{16})$/.test(req.body.expId)) {
	var expID = req.body.expId;
	var pid=expID.slice(0,4);
	var expUnique=expID.slice(4,20);
        var expId = {expId : expID };
	if(fs.existsSync('experimentData/'+pid+'/'+expUnique+'/')) {
	    var proteinInfo;
	    try{
		data = fs.readFileSync(path.join('/home/completePocketsNRSet',pid + 'A' + '.csv'), {encoding : 'utf8'});
		var options = {
                    delimiter : " ",
                    wrap : true,
                    headers :"pocket area residues"};
		pockets = csvToJson.toObject(data,options);
		console.log(pockets);
		proteinInfo = {pid:pid,expId, pockets};
	    } catch (err) {
		proteinInfo = {pid:pid,expId};
		console.log('surf report doesnt exist');
	    }
	    res.send(JSON.stringify(proteinInfo));
	} else {
	    // this needs handling
	    res.send('invalid')
	}
    } else {
	// see above, since frontend cannot send anything invalid
	res.send('backend experiment id check failed!');
    }
});	

router.post('/getPidsWithSurfReports',(req,res) => {
	
	var data = fs.readFileSync('./A_ChainPDBs.txt').toString().split('\n');
	console.log(data);
	res.send(data);
});
	



module.exports = router;
