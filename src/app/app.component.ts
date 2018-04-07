import { Component } from '@angular/core';
import {Router} from '@angular/router';
import {ServerService} from './server.service';
import {UserInputService} from './userinput.service';
@Component({
  selector: 'app-root',
  templateUrl: './app.component.html',
  styleUrls: ['./app.component.css']
})
export class AppComponent {

  constructor(private serverService: ServerService,
  private router: Router,private userInputService: UserInputService ){
    this.userInputService.pidEntered.subscribe(
       (pid : string) => this.submitPidToServer(pid));
   }

  title = 'Protein Cavity Explorer';
  experiments = [
    {
    pid: '32',
    name: 'david'}
  ];


  onSubmit(){
    this.serverService.createExperiment({pid: 'lake'})
    .subscribe((response) => {
      console.log(response);
      console.log(response['pid'])},
    (error) => {console.log(error)});
  }

  submitPidToServer(pid : string){
    console.log("submitPidToServer " + pid);
    this.router.navigate(['experiment']);
  }


  submitExperimentId(experimentId : string){
    console.log(experimentId);

  }

}
