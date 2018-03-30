import { Component } from '@angular/core';
import {ServerService} from './server.service';
@Component({
  selector: 'app-root',
  templateUrl: './app.component.html',
  styleUrls: ['./app.component.css']
})
export class AppComponent {
  constructor(private ServerService: ServerService){}
  title = 'Protein Cavity Explorer';
  experiments = [
    {
    pid: '32',
    name: 'david'}
  ];


  onSubmit(){
    this.ServerService.createExperiment({pid: 'lake'})
    .subscribe((response) => {
      console.log(response);
      console.log(response['pid'])},
    (error) => {console.log(error)});
  }
}
